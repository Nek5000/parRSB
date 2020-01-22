#include <stdlib.h>
#include <math.h>

#include "gencon-impl.h"
#include "exa-memory.h"
#include "exasort.h"

int mergeSegments(exaHandle h,Mesh mesh,int i,exaScalar tolSquared)
{
  exaComm c  =exaGetComm(h);

  Point points=exaArrayGetPointer(mesh->elements);
  exaInt nPoints=exaArrayGetSize(mesh->elements);

  exaInt bin=1;
  if(nPoints==0) bin=0;
  exaComm nonZeroRanks;
  exaCommDup(&nonZeroRanks,c);
  exaCommSplit(&nonZeroRanks,bin);
  exaInt rank=exaCommRank(nonZeroRanks);
  exaInt size=exaCommSize(nonZeroRanks);

  if(bin==1){
    exaArray arr;
    exaArrayInit(&arr,struct Point_private,1);

    Point lastp; exaMalloc(1,&lastp);
    memcpy(lastp,&points[nPoints-1],sizeof(*lastp));
    lastp->proc=(rank+1)%size;
    exaArrayAppend(arr,lastp);
    exaFree(lastp);

    exaArrayTransfer(arr,offsetof(struct Point_private,proc),1,
      nonZeroRanks);

    exaInt n=exaArrayGetSize(arr);
    assert(n==1);
    lastp=exaArrayGetPointer(arr);

    if(rank>0){
      exaScalar d=sqrDiff(lastp->x[i],points->x[i]);
      exaScalar dx=min(lastp->dx,points->dx)*tolSquared;
      if(d>dx) points->ifSegment=1;
    }

    exaArrayFree(arr);
  }

  exaDestroy(nonZeroRanks);
  return 0;
}

int findLocalSegments(Mesh mesh,int i,exaScalar tolSquared){
  Point      ptr=exaArrayGetPointer(mesh->elements);
  exaInt nPoints=exaArrayGetSize(mesh->elements);

  exaInt j;
  for(j=0;j<nPoints-1;j++){
    exaScalar d=sqrDiff(ptr[j].x[i],ptr[j+1].x[i]);
    exaScalar dx=min(ptr[j].dx,ptr[j+1].dx)*tolSquared;
    if(d>dx) ptr[j+1].ifSegment=1;
  }
}

int findSegments(exaHandle h,Mesh mesh,exaScalar tol){
  int nDim=mesh->nDim;
  int nVertex=mesh->nVertex;
  exaInt rank=exaRank(h);
  exaScalar tolSquared=tol*tol;

  exaSort(mesh->elements,
    exaScalar_t,offsetof(struct Point_private,x[0]),
    exaSortAlgoBinSort,0,exaGetComm(h));

  Point   points=exaArrayGetPointer(mesh->elements);
  exaInt nPoints=exaArrayGetSize(mesh->elements);
  exaInt i;

  exaSortArray3(mesh->elements,
    exaScalar_t,offsetof(struct Point_private,x[0]),
    exaScalar_t,offsetof(struct Point_private,x[1]),
    exaScalar_t,offsetof(struct Point_private,x[2]));

  exaLong out[2][1],buff[2][1],in[1];
  in[0]=nPoints;
  exaCommScan(exaGetComm(h),out,in,buff,1,exaLong_t,exaAddOp);
  exaLong start=out[0][0];

  for(i=0;i<nPoints;i++) points[i].globalId=start+i;

  exaSort(mesh->elements,
    exaLong_t,offsetof(struct Point_private,globalId),
    exaSortAlgoBinSort,1,exaGetComm(h));

  points=exaArrayGetPointer(mesh->elements);
  nPoints=exaArrayGetSize(mesh->elements);
  exaDebug(h,"rank=%d load=%d\n",rank,nPoints);

  for(i=0;i<nPoints;i++) points[i].ifSegment=0,points[i].proc=rank;

  findLocalSegments(mesh,0,tolSquared);
  mergeSegments  (h,mesh,0,tolSquared);

  int ipass;
  for(ipass=0; ipass<nDim; ipass++)
    for(i=0;i<nDim;i++){
      findLocalSegments(mesh,i,tolSquared);
      mergeSegments  (h,mesh,i,tolSquared);
    }

  return 0;
}

#undef sqrDiff
#undef min
