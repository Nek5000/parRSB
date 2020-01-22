#include "gencon-impl.h"

int neighborMap[GC_MAX_VERTICES][GC_MAX_NEIGHBORS]={
  {1,2,4},
  {0,3,5},
  {0,3,6},
  {1,2,7},
  {0,5,6},
  {1,4,7},
  {2,4,7},
  {3,5,6}
};

int findMinNeighborDistance(exaHandle h,Mesh mesh){
  Element p=exaArrayGetPointer(mesh->elements),e;
  e=p+mesh->nelt;

  int nDim=mesh->nDim;
  int nVertex=mesh->nVertex;
  int nNeighbors=mesh->nNeighbors;

  int j,k,neighbor;
  exaScalar d;
  for(; p!=e; p++){
    for(j=0; j<nVertex; j++){
      p->vertex[j].dx=exaScalarMAX;
      for(k=0; k<nNeighbors; k++){
        neighbor=neighborMap[j][k];
        if(nDim==3) d=distance3D(p->vertex[j],p->vertex[neighbor]);
        else d=distance2D(p->vertex[j],p->vertex[neighbor]);
        p->vertex[j].dx=min(p->vertex[j].dx,d);
      }
    }
  }
}
