#include <genmap-impl.h>
#include <exa-memory.h>

#define min(a,b) ((b)<(a)?(b):(a))

typedef struct{
  GenmapLong sequenceId;
  GenmapLong elementId;
  GenmapLong neighbors[8];
  int nNeighbors;
  GenmapLong vertexId;
  int workProc;
} vertex;

typedef struct{
  GenmapULong elementId;
} element;

int GenmapFindNeighbors(GenmapHandle h,GenmapComm c,GenmapLong **eIds_,
    GenmapInt **offsets_)
{
  GenmapInt lelt=GenmapGetNLocalElements(h);
  GenmapInt nv  =GenmapGetNVertices(h);

  GenmapScan(h,c);
  GenmapLong elementId =GenmapGetLocalStartIndex(h)+1;
  GenmapLong sequenceId=elementId*nv;

#if defined(GENMAP_DEBUG)
  printf("nid=%d startId=%lld\n",GenmapCommRank(c),elementId);
#endif

  size_t size=lelt*nv;
  exaArray vertices; exaArrayInit(&vertices,vertex,size);

  GenmapElements elems=GenmapGetElements(h);
  GenmapInt i; int j;
  for(i=0;i<lelt;i++){
    for(j=0;j<nv;j++){
      vertex t={
        .elementId =elementId,
        .sequenceId=sequenceId,
        .vertexId  =elems[i].vertices[j],
        .workProc  =elems[i].vertices[j]%GenmapCommSize(c)
      };
#if defined(GENMAP_DEBUG)
      printf("vid=%lld workProc=%d\n",t.vertexId,t.workProc);
#endif
      exaArrayAppend(vertices,(void*)&t);
      sequenceId++;
    }
    elementId++;
  }

  exaComm comm; exaCommCreate(&comm,c->gsc.c);
  exaArrayTransfer(vertices,offsetof(vertex,workProc),1,comm);

  assert(sizeof(GenmapLong)==sizeof(exaLong));
  assert(sizeof(GenmapInt)==sizeof(exaInt));

  exaSortArray(vertices,exaLong_t,offsetof(vertex,vertexId));

  size=exaArrayGetSize(vertices);
  vertex* vPtr=exaArrayGetPointer(vertices);
#if defined(GENMAP_DEBUG)
  printf("nid=%d lelt=%lu ",exaCommRank(comm),size);
  for(i=0;i<size;i++)
    printf("%lld ",vPtr[i].vertexId);
  printf("\n");
#endif

  //FIXME: Assumes quads or hexes
  exaInt s=0,e;
  while(s<size){
    for(e=s+1; e<size && vPtr[s].vertexId==vPtr[e].vertexId; e++);
    int nNeighbors=min(e,size)-s;
    for(i=s;i<min(e,size);i++){
      for(j=0;j<nNeighbors;j++)
        vPtr[i].neighbors[j]=vPtr[s+j].elementId;
      vPtr[i].nNeighbors=nNeighbors;
    }
    s=e;
  }

  exaArrayTransfer(vertices,offsetof(vertex,workProc),0,comm);

  exaSortArray(vertices,exaLong_t,offsetof(vertex,sequenceId));
  vPtr=exaArrayGetPointer(vertices);
  size=exaArrayGetSize(vertices); assert(size==lelt*nv); // sanity check

#if defined(GENMAP_DEBUG)
  for(i=0;i<size;i++)
    printf("vid=%lld neighbors=%d\n",vPtr[i].vertexId,vPtr[i].nNeighbors);
#endif

  exaArray nbrs;
  exaArrayInit(&nbrs,element,GENMAP_MAX_VERTICES*GENMAP_MAX_NEIGHBORS);

  exaMalloc(lelt*27,eIds_   ); exaLong *eIds   =*eIds_;
  exaMalloc(lelt+ 1,offsets_); exaInt  *offsets=*offsets_;

  exaInt cnt=0; int k; offsets[lelt]=0;
  for(i=0;i<lelt;i++){
    exaArraySetSize(nbrs,0);

    element e; GenmapLong curId;
    curId=e.elementId=vPtr[i*nv].elementId;
    for(j=0;j<nv;j++){
      vertex v=vPtr[i*nv+j];
      for(k=0;k<v.nNeighbors;k++)
        if((e.elementId=v.neighbors[k])!=curId)
          exaArrayAppend(nbrs,(void*)&e);
    }

    exaSortArray(nbrs,exaULong_t,offsetof(element,elementId));
    element *ePtr=exaArrayGetPointer(nbrs);
    size=exaArrayGetSize(nbrs);

#if defined(GENMAP_DEBUG)
    printf("neighbors: ");
    for(k=0;k<size;k++)
      printf("%lld ",ePtr[k].elementId);
    printf("\n");
#endif

    offsets[i]=0;
    for(eIds[cnt++]=curId,j=0; j<size; j++)
      if(eIds[cnt-1]!=-ePtr[j].elementId)
        eIds[cnt++]=-ePtr[j].elementId,eIds[i]+=1.0;
    offsets[lelt]+=offsets[i];
  }

  exaDestroy(nbrs);
  exaDestroy(vertices);
  exaCommDestroy(comm);
}

int GenmapInitLaplacian(GenmapHandle h,GenmapComm c,GenmapVector weights)
{
  GenmapLong *eIds; GenmapInt *offsets;
  GenmapFindNeighbors(h,c,&eIds,&offsets);

  GenmapInt lelt=GenmapGetNLocalElements(h);
  GenmapInt nv  =GenmapGetNVertices(h);

  GenmapInt i;
  for(i=0;i<lelt;i++)
    weights->data[i]=offsets[i];

#if defined(GENMAP_DEBUG)
  printf("weights: ");
  for(i=0;i<lelt;i++){
    printf(" (%lld,%lf)",vPtr[i*nv].elementId,weights->data[i]);
  }
  printf("\n");
  int cnt1=0;
  for(i=0;i<lelt;i++){
    printf("%lld:",vPtr[i*nv].elementId);
    for(int j=0;j<weights->data[i]+1;j++)
      printf(" %lld",eIds[cnt1++]);
    printf("\n");
  }
  printf("\n");
#endif

  c->gsh=gs_setup(eIds,offsets[lelt],&c->gsc,0,gs_crystal_router,0);
  GenmapMalloc(offsets[lelt],&c->laplacianBuf);

  exaFree(eIds);
  exaFree(offsets);

  return 0;
}

int GenmapLaplacian(GenmapHandle h,GenmapComm c,GenmapVector u,
  GenmapVector weights,GenmapVector v)
{
  GenmapInt lelt=GenmapGetNLocalElements(h);
  GenmapInt nv  =GenmapGetNVertices(h);

  assert(u->size==v->size);
  assert(u->size==lelt   );

  GenmapInt i,cnt; int j;
  for(i=0,cnt=0; i<lelt; cnt+=weights->data[i],i++)
    c->laplacianBuf[cnt++]=u->data[i];

  gs(c->laplacianBuf,genmap_gs_scalar,gs_add,0,c->gsh,
    &c->buf);

  for(cnt=0,i=0; i<lelt; i++){
    int nbrs=weights->data[i];
    v->data[i]=c->laplacianBuf[cnt++]*nbrs;
    for(j=0; j<nbrs; j++)
      v->data[i]-=c->laplacianBuf[cnt++];
  }

  return 0;
}

#undef min
