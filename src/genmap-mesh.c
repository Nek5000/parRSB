#include <genmap-impl.h>
#include <gencon-impl.h>

#include <exa-memory.h>

int MeshInit(Mesh *m_,int nel,int nDim){
  exaMalloc(1,m_);
  Mesh m=*m_;

  m->nelt=nel;
  m->nDim=nDim;
  m->nNeighbors=nDim;
  m->nVertex=(nDim==2)?4:8;

  exaArrayInit(&m->elements,struct Element_private ,nel);
  exaArrayInit(&m->boundary,struct Boundary_private,  0);

  return 0;
}

Element MeshGetElements(Mesh m){
  return exaArrayGetPointer(m->elements);
}

int MeshFree(Mesh m){
  exaArrayFree(m->elements);
  exaArrayFree(m->boundary);

  exaFree(m);

  return 0;
}
