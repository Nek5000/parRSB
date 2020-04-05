#include "exa-impl.h"
#include "exa-types.h"
#include "exasort.h"
#include "parRSB.h"

void getAxisLength(exaScalar **length,exaArray eList,exaComm comm,
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

  exaCommGop(comm,min,3,exaScalar_t,exaMinOp);
  exaCommGop(comm,max,3,exaScalar_t,exaMaxOp);

  exaMalloc(MAXDIM,length);

  for(i=0;i<ndim;i++)
    (*length)[i]=max[i]-min[i];
}

int parRCB(exaComm comm,exaArray eList,int ndim){
  exaInt rank=exaCommRank(comm);
  exaInt size=exaCommSize(comm);

  while(size>1){
    exaScalar *length=NULL;
    getAxisLength(&length,eList,comm,ndim);

    int axis1=0,axis3=0,d;
    for(d=1;d<ndim;d++){
      if(length[axis1]<length[d]) axis1=d;
      if(length[axis3]>length[d]) axis3=d;
    }

    exaUInt offsets[3]={offsetof(elm_rcb,coord[0]),
      offsetof(elm_rcb,coord[1]),offsetof(elm_rcb,coord[2])};

    exaSort2(eList,exaScalar_t,offsets[axis1],exaScalar_t,
      offsets[axis3],exaSortAlgoBinSort,1,comm);

    int p=(size+1)/2;
    int bin=(rank>=p);
    exaCommSplit(&comm,bin,rank);

    exaFree(length);

    rank=exaCommRank(comm);
    size=exaCommSize(comm);
  }

  return 0;
}
