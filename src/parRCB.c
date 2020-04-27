#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "exa-impl.h"
#include "exasort.h"
#include "parRSB.h"

void fparRCB_partMesh(int *part,double *vtx,int *nel,int *nv,
  int *options,int *comm,int *err)
{
  *err = 1;

  exaCommExternal c;
  c = MPI_Comm_f2c(*comm);
  *err=parRCB_partMesh(part,vtx,*nel,*nv,options,c);
}

int parRCB_partMesh(int *part,double *vtx,int nel,int nv,
  int *options,MPI_Comm comm)
{
  exaHandle h;
  exaInit(&h,comm,"/host");

  int rank=exaRank(h),size=exaSize(h);

  /* load balance input data */
  exaLong nelg=0;
  exaLong nell=nel;
  exaAllReduce(h,&nelg,&nell,1,exaLong_t,exaAddOp);

  exaLong nelg_start,buf0;
  exaScan(h,&nelg_start,&nell,&buf0,1,exaLong_t,exaAddOp);

  exaArray eList; exaArrayInit(&eList,elm_rcb,nel);
  elm_rcb *data=exaArrayGetPointer(eList);

  int ndim=(nv==8)?3:2;

  int e, n;
  for(e=0;e<nel;++e){
    data[e].id=nelg_start+(e+1);
    data[e].orig=rank;
    for(int n=0;n<ndim;n++)
      data[e].coord[n]=vtx[e*ndim+n];
  }

  exaArraySetSize(eList,nel);

  exaLoadBalance(eList,exaGetComm(h));

  nel=exaArrayGetSize(eList);
  exaComm commRcb;
  exaCommSplit(exaGetComm(h),nel>0,rank,&commRcb);

  if(nel>0){
    double time0 = comm_time();

    parRCB(commRcb,eList,ndim);

    if(exaRank(h)==0)
      printf("\nparRCB finished in %lfs\n",comm_time()-time0);

    fflush(stdout);
  }

  /* restore original input */
  exaArrayTransfer(eList,offsetof(elm_rcb,orig),1,exaGetComm(h));
  exaSortArray(eList,exaLong_t,offsetof(elm_rcb,id));

  exaBarrier(h);

  data=exaArrayGetPointer(eList);
  nel=exaArrayGetSize(eList);

  for(int e=0;e<nel;e++) {
    part[e]=data[e].orig;
  }

  exaArrayFree(eList);

  exaFinalize(h);

  return 0;
}
