#include <genmap-impl.h>
#include <exa-memory.h>

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
  size_t size   =lelt*nv;
  exaArray vertices; exaArrayInit(&vertices,vertex,size);

  int np=GenmapCommSize(c);

  GenmapScan(h,c);
  GenmapLong elementId =GenmapGetLocalStartIndex(h);
  GenmapLong sequenceId=elementId*nv;

  GenmapElements elems=GenmapGetElements(h);

  GenmapInt i;
  int j;
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
    for(e=s+1; e<size && vPtr[e].vertexId==vPtr[s].vertexId; e++);
    if(e<size){
      exaInt nNeighbors=e-s;
      // vertex can't be shared by more than 8 elements
      assert(nNeighbors<=GENMAP_MAX_NEIGHBORS);
      for(i=s;i<e;i++){
        for(j=0;j<nNeighbors;j++){
          vPtr[i].neighbors[j]=vPtr[s+j].elementId;
        }
        vPtr[i].nNeighbors=nNeighbors;
      }
    }
    s=e;
  }

  exaArrayTransfer(vertices,offsetof(vertex,workProc),0,comm);
  exaSortArray(vertices,exaLong_t,offsetof(vertex,sequenceId));

  vPtr=exaArrayGetPointer(vertices);
  size=exaArrayGetSize(vertices); assert(size==lelt*nv); // sanity check

  exaArray neighbors; exaArrayInit(&neighbors,element,\
    GENMAP_MAX_VERTICES*GENMAP_MAX_NEIGHBORS);

  exaLong *elementIds; exaMalloc(lelt*nv,&elementIds);

  exaInt count=0;
  int k;
  for(i=0;i<lelt;i++){
    weights->data[i]=0.0;

    exaArraySetSize(neighbors,0);
    for(j=0;j<nv;j++){
      vertex v=vPtr[i*nv+j];
      element e;
      for(k=0;k<v.nNeighbors;k++)
        e.elementId=v.neighbors[k],exaArrayAppend(neighbors,(void*)&e);
    }
    size=exaArrayGetSize(neighbors);

    exaSortArray(neighbors,exaLong_t,offsetof(element,elementId));

    element *ePtr=exaArrayGetPointer(neighbors);
    elementIds[count++]=ePtr[0].elementId,weights->data[i]+=1.0;
    for(j=1;j<size;j++)
      if(elementIds[count-1]!=ePtr[j].elementId)
        elementIds[count++]=ePtr[j].elementId,weights->data[i]+=1.0;
  }
  exaFree(neighbors);
  exaFree(vertices);

  c->laplacianHandle=gs_setup(elementIds,count,&c->gsComm,0,
    gs_crystal_router,0);

  exaFree(elementIds);
  exaCommDestroy(comm);
}
