#include "exa-impl.h"
#include "exa-types.h"
#include "exasort.h"
#include "parRSB.h"

void getAxisLength(exaScalar *length,exaArray eList,exaComm comm,
  int ndim)
{
  exaScalar min[MAXDIM],max[MAXDIM];
  exaInt i;
  for(i=0;i<ndim;i++) min[i]=exaScalarMAX,max[i]=exaScalarMIN;

  exaInt nel=exaArrayGetSize(eList);
  elm_rcb* elems=exaArrayGetPointer(eList);
  for(i=0;i<nel;i++){
    if(elems[i].coord[0]<min[0]) min[0]=elems[i].coord[0];
    if(elems[i].coord[1]<min[1]) min[1]=elems[i].coord[1];
    if(ndim==3)
      if(elems[i].coord[2]<min[2]) min[2]=elems[i].coord[2];

    if(elems[i].coord[0]>max[0]) max[0]=elems[i].coord[0];
    if(elems[i].coord[1]>max[1]) max[1]=elems[i].coord[1];
    if(ndim==3)
      if(elems[i].coord[2]>max[2]) max[2]=elems[i].coord[2];
  }

  exaGop(comm,min,3,exaScalar_t,exaMinOp);
  exaGop(comm,max,3,exaScalar_t,exaMaxOp);

  exaMalloc(2*MAXDIM,&length);
  for(i=0;i<ndim;i++){
    length[2*i+0]=max[i]-min[i];
    length[2*i+1]=0.5*(max[i]+min[i]);
  }
}

int parRCB(exaComm comm,exaArray eList,int ndim){
  exaInt rank=exaCommRank(comm);
  exaInt size=exaCommSize(comm);

  while(size>1){
    exaScalar *length;
    getAxisLength(length,eList,comm,ndim);

    int axis=0,d;
    for(d=1;d<ndim;d++)
      if(length[2*axis]<length[2*d]) axis=d;
    printf("size: %d x: %lf %f,y: %lf %lf,axis=",
      size,length[0],length[1],length[2],length[3],axis);

    int nel=exaArrayGetSize(eList);
    elm_rcb *elems=exaArrayGetPointer(eList);

    exaInt low=0,high=0;
    for(d=0;d<nel;d++)
      if(elems[d].coord[axis]<length[2*axis+1]) low++;
      else high++;

    exaLong in[2],out[2][2],buf[2][2];
    in[0]=low,in[1]=high;
    exaCommScan(comm,out,in,buf,2,exaLong_t,exaAddOp);
    exaLong lowStart=out[0][0],nGlobalLow=out[1][0];
    exaLong highStart=out[0][1],nGlobalHigh=out[1][1];

    printf("low=%d lowStarrt=%d high=%d highStart=%d\n",
      low,lowStart,high,highStart);

    switch(axis){
      case 0:
        exaSortArray(eList,exaScalar_t,
          offsetof(elm_rcb,coord[0]));
        break;
      case 1:
        exaSortArray(eList,exaScalar_t,
          offsetof(elm_rcb,coord[1]));
        break;
      case 2:
        exaSortArray(eList,exaScalar_t,
          offsetof(elm_rcb,coord[2]));
        break;
      default:
        break;
    }

    int p=size*nGlobalLow/(nGlobalLow+nGlobalHigh);
    for(d=0;d<low;d++){
      elems[d].proc=(lowStart+d)*p/nGlobalLow;
    }
    for(d=0;d<high;d++){
      elems[d+low].proc=p+(highStart+d)*(size-p)/nGlobalHigh;
    }
    exaArrayTransfer(eList,offsetof(elm_rcb,proc),0,comm);

    int bin=(rank>=p);
    exaCommSplit(&comm,bin,rank);

    exaFree(length);

    parRCB(comm,eList,ndim);

    rank=exaCommRank(comm);
    size=exaCommSize(comm);
  }

  return 0;
}
