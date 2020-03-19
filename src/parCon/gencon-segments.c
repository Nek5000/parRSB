#include <stdlib.h>
#include <math.h>

#include "gencon-impl.h"
#include "exa-memory.h"
#include "exasort.h"

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

int sortInsideSegments(Mesh mesh,int dim){
  exaInt nPoints=exaArrayGetSize(mesh->elements);
  Point end=(Point)exaArrayGetPointer(mesh->elements)+nPoints-1;
  Point lstart=exaArrayGetPointer(mesh->elements),lend;

  buffer buf; buffer_init(&buf,1024);

  while(lstart!=end){
    lend=lstart+1;
    while(lend!=end && lend->ifSegment==0) lend++;

    exaInt nLocal=lend-lstart;
    if(nLocal>1) {
      switch(dim){
        case 0:
          sarray_sort  (struct Point_private,lstart,nLocal,
            x[0],3,&buf);
          break;
        case 1:
          sarray_sort_2(struct Point_private,lstart,nLocal,
            x[0],3,x[1],3,&buf);
          break;
        case 2:
          sarray_sort_3(struct Point_private,lstart,nLocal,
            x[0],3,x[1],3,x[2],3,&buf);
          break;
        default:
          break;
      }
    }
    lstart=lend;
  }

  buffer_free(&buf);
}

int mergeSegments(exaHandle h,Mesh mesh,int i,
  exaScalar tolSquared)
{
  exaComm c=exaGetComm(h);

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

int findSegments(exaHandle h,Mesh mesh,exaScalar tol){
  int nDim=mesh->nDim;
  int nVertex=mesh->nVertex;
  exaInt rank=exaRank(h);
  exaScalar tolSquared=tol*tol;

  exaSort(mesh->elements,
    exaScalar_t,offsetof(struct Point_private,x[0]),
    exaSortAlgoBinSort,0,exaGetComm(h));

  Point points=exaArrayGetPointer(mesh->elements);
  exaInt nPoints=exaArrayGetSize(mesh->elements);

  exaLong out[2][1],buff[2][1],in[1];
  in[0]=nPoints;
  exaCommScan(exaGetComm(h),out,in,buff,1,exaLong_t,exaAddOp);
  exaLong start=out[0][0];
  exaLong nGlobalPoints=out[1][0];

  exaInt i;
  for(i=0;i<nPoints;i++){
    points[i].ifSegment=0;
    points[i].proc=rank;
  }

  int ipass,dim;
  for(ipass=0;ipass<nDim;ipass++){
    for(dim=0;dim<nDim;dim++){
      if(exaGetDebug(h)>0){
        exaInt count=0;
        for(i=0;i<nPoints;i++){
          if(points[i].ifSegment>0){
            count++;
          }
        }

        exaLong in[1];
        in[0]=count;
        exaCommGop(exaGetComm(h),in,1,exaLong_t,exaAddOp);
        if(exaRank(h)==0)
          printf("locglob: %d %lld %lld\n",dim+1,in[0]+1,\
            nGlobalPoints);
      }

      findLocalSegments(mesh,dim,tolSquared);
      mergeSegments(h,mesh,dim,tolSquared);
      if(dim<2) sortInsideSegments(mesh,dim+1);
    }
  }

  //for(i=0;i<nPoints;i++){
  //  printf("rank=%d i=%02d %d %lf %lf %lf\n",rank,
  //    i+1,points[i].ifSegment,points[i].x[0],points[i].x[1],
  //    points[i].x[2]);
  //}

  return 0;
}

#undef sqrDiff
#undef min
