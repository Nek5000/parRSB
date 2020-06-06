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
  GenmapLong elementId;
} element;

int GenmapInitLaplacian(GenmapHandle h,GenmapComm c,GenmapVector weights)
{
  GenmapInt lelt=GenmapGetNLocalElements(h);
  GenmapInt nv  =GenmapGetNVertices(h);

  int np=GenmapCommSize(c);
  GenmapScan(h,c);
  GenmapLong elementId =GenmapGetLocalStartIndex(h);
  GenmapLong sequenceId=elementId*nv;

  size_t size=lelt*nv;
  exaArray vertices; exaArrayInit(&vertices,vertex,size);

  GenmapElements elems =GenmapGetElements(h);
  GenmapInt i; int j;
  for(i=0;i<lelt;i++){
    for(j=0;j<nv;j++){
      vertex t={
        .elementId =elementId,
        .sequenceId=sequenceId,
        .vertexId  =elems[i].vertices[j],
        .workProc  =elems[i].vertices[j]%np
      };
      exaArrayAppend(vertices,(void*)&t);
      sequenceId++;
    }
    elementId++;
  }

  exaComm comm; exaCommCreate(&comm,c->gsComm.c);
  exaArrayTransfer(vertices,offsetof(vertex,workProc),1,comm);

  assert(sizeof(GenmapLong)==sizeof(exaLong));
  exaSortArray(vertices,exaLong_t,offsetof(vertex,vertexId));

  size=exaArrayGetSize(vertices);
  vertex* vPtr=exaArrayGetPointer(vertices);

  //TODO: Assumes quads or hexes
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

  exaArray neighbors;
  exaArrayInit(&neighbors,element,\
    GENMAP_MAX_VERTICES*GENMAP_MAX_NEIGHBORS);

  exaLong *elementIds; exaMalloc(lelt*27,&elementIds);

  exaInt count=0; int k;
  for(i=0;i<lelt;i++){
    weights->data[i]=0.0;

    exaArraySetSize(neighbors,0);
    for(j=0;j<nv;j++){
      vertex v=vPtr[i*nv+j];
      element e;
      for(k=0;k<v.nNeighbors;k++){
        e.elementId=v.neighbors[k];
        exaArrayAppend(neighbors,(void*)&e);
      }
    }

    exaSortArray(neighbors,exaLong_t,offsetof(element,elementId));
    size=exaArrayGetSize(neighbors); assert(size>0);

    element *ePtr=exaArrayGetPointer(neighbors);
    elementIds[count]=ePtr[0].elementId,count++;
    for(j=1;j<size;j++)
      if(elementIds[count-1]!=ePtr[j].elementId)
        elementIds[count]=ePtr[j].elementId,weights->data[i]+=1.0,count++;
  }
  exaDestroy(neighbors);
  exaDestroy(vertices );

  c->laplacianHandle=gs_setup(elementIds,count,&c->gsComm,0,
    gs_crystal_router,0);

  exaFree(elementIds);
  exaCommDestroy(comm);
}

#undef min
