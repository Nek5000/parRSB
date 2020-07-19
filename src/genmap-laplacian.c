#include <genmap-impl.h>
#include <exa-memory.h>

#define min(a,b) ((b)<(a)?(b):(a))

typedef struct{
  GenmapULong sequenceId;
  GenmapLong neighbors[8];
  int nNeighbors;
  GenmapULong elementId;
  GenmapULong vertexId;
  uint workProc;
} vertex;

typedef struct{
  GenmapULong elementId;
} element;

int GenmapFindNeighbors(GenmapHandle h,GenmapComm c,GenmapLong **eIds_,
    GenmapInt **neighbors_)
{
  struct comm cc=c->gsc;

  sint lelt=GenmapGetNLocalElements(h);
  sint nv  =GenmapGetNVertices(h);

  GenmapScan(h,c);
  ulong elem_id   =GenmapGetLocalStartIndex(h)+1;
  ulong sequenceId=elem_id*nv;

  size_t size=lelt*nv;
  struct array vertices; array_init(vertex,&vertices,size);

  GenmapElements elems=GenmapGetElements(h);
  sint i,j;
  for(i=0;i<lelt;i++){
    for(j=0;j<nv;j++){
      vertex t={
        .elementId =elem_id,
        .sequenceId=sequenceId,
        .vertexId  =elems[i].vertices[j],
        .workProc  =elems[i].vertices[j]%cc.np
      };
      array_cat(vertex,&vertices,&t,1);
      sequenceId++;
    }
    elem_id++;
  }

  struct crystal cr; crystal_init(&cr,&cc);
  sarray_transfer(vertex,&vertices,workProc,1,&cr);
  size=vertices.n; vertex* vPtr=vertices.ptr;

  buffer buf; buffer_init(&buf,1024);
  sarray_sort(vertex,vPtr,size,vertexId,1,&buf);

  struct array a; array_init(csr_entry,&a,10);

  //FIXME: Assumes quads or hexes
  exaInt s=0,e; csr_entry t;
  while(s<size){
    for(e=s+1; e<size && vPtr[s].vertexId==vPtr[e].vertexId; e++);
    int nNeighbors=min(e,size)-s;
    for(i=s;i<min(e,size);i++){
      // get rid of these
      for(j=0;j<nNeighbors;j++)
        vPtr[i].neighbors[j]=vPtr[s+j].elementId;
      vPtr[i].nNeighbors=nNeighbors;

      t.r=vPtr[i].elementId; t.proc=vPtr[i].workProc;
      for(j=0;j<nNeighbors;j++){
        t.c=vPtr[s+j].elementId;
        array_cat(csr_entry,&a,&t,1);
      }
    }
    s=e;
  }

  sarray_transfer(vertex,&vertices,workProc,0,&cr);
  vPtr=vertices.ptr; size=vertices.n; assert(size==lelt*nv);//sanity-check

  sarray_transfer(csr_entry,&a,proc,1,&cr);
  crystal_free(&cr);

  sarray_sort(vertex,vertices.ptr,vertices.n,sequenceId,1,&buf);
  sarray_sort_2(csr_entry,a.ptr,a.n,r,1,c,1,&buf);

  buffer_free(&buf);

  exaArray nbrs;
  exaArrayInit(&nbrs,element,GENMAP_MAX_VERTICES*GENMAP_MAX_NEIGHBORS);

  exaMalloc(lelt*27,eIds_   ); exaLong *eIds   =*eIds_;
  exaMalloc(lelt+ 1,neighbors_); exaInt  *neighbors=*neighbors_;

  exaInt cnt=0; int k; neighbors[lelt]=0;
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

    neighbors[i]=0;
    for(eIds[cnt++]=curId,j=0; j<size; j++)
      if(eIds[cnt-1]!=-ePtr[j].elementId)
        eIds[cnt++]=-ePtr[j].elementId,neighbors[i]+=1;
    neighbors[lelt]+=neighbors[i]+1;
  }

  array_free(&vertices);
  exaDestroy(nbrs);
}

int GenmapInitLaplacian(GenmapHandle h,GenmapComm c,GenmapVector weights)
{
  GenmapLong *eIds; GenmapInt *neighbors;
  GenmapFindNeighbors(h,c,&eIds,&neighbors);

  GenmapInt lelt=GenmapGetNLocalElements(h);
  GenmapInt nv  =GenmapGetNVertices(h);

  GenmapInt i;
  for(i=0;i<lelt;i++)
    weights->data[i]=neighbors[i];

  c->gsh=gs_setup(eIds,neighbors[lelt],&c->gsc,0,gs_crystal_router,0);
  GenmapMalloc(neighbors[lelt],&c->laplacianBuf);

  exaFree(eIds);
  exaFree(neighbors);

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
