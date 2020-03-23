#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "exa-impl.h"
#include "exasort.h"
#include "parRSB.h"

void fparRCB_partMesh(int *part,double *vtx,int *nel,int *ndim,
  int *options,int *comm,int *err)
{
  *err = 1;

  exaCommExternal c;
  c = MPI_Comm_f2c(*comm);
  *err=parRCB_partMesh(part,vtx,*nel,*ndim,options,c);
}

int parRCB_partMesh(int *part,double *vtx,int nel,int ndim,
  int *options,MPI_Comm comm)
{
  exaHandle h;
  exaInit(&h,comm,"/host");

  int rank=exaRank(h),size=exaSize(h);

  /* load balance input data */
  exaLong nelg=0;
  exaLong nell=nel;
  exaAllReduce(h,&nelg,&nell,1,exaLong_t,exaAddOp);
  exaLong nstar=nelg/size;
  if(nstar==0) nstar=1;
  if(rank==0)
    printf("nelg: %lld nell: %lld nstar: %lld\n",nelg,nell,nstar);

  exaLong nelg_start,buf0;
  exaScan(h,&nelg_start,&nell,&buf0,1,exaLong_t,exaAddOp);

  exaArray eList; exaArrayInit(&eList,elm_rcb,nel);
  elm_rcb *data=exaArrayGetPointer(eList);

  int e, n;
  for(e=0;e<nel;++e){
    data[e].id=nelg_start+(e+1);
    const exaLong eg=data[e].id;

    data[e].proc=(int) ((eg-1)/nstar);
    data[e].orig=rank;
    if(eg>size*nstar) data[e].proc= (eg%size)-1; 
    for(int n=0;n<ndim;n++)
      data[e].coord[n]=vtx[e*ndim+n];
  }
  exaArraySetSize(eList,nel);

  buffer buf; buffer_init(&buf,1024);

  exaArrayTransfer(eList,offsetof(elm_rcb,proc),1,exaGetComm(h));
  exaSortArray(eList,exaLong_t,offsetof(elm_rcb,id));

  exaComm commRcb; exaCommDup(&commRcb,exaGetComm(h));
  nel=exaArrayGetSize(eList);
  exaCommSplit(&commRcb,nel>0,rank);

  if(nel>0){
    double time0 = comm_time();
    if(options!=NULL) exaSetDebug(h,options[0]);
    if(exaRank(h)==0)
      exaDebug(h,"\nfinished in %lfs\n",comm_time()-time0);

    parRCB(commRcb,eList,ndim);

    fflush(stdout);
  }

  /* restore original input */
  exaArrayTransfer(eList,offsetof(elm_rcb,orig),1,exaGetComm(h));
  exaSortArray(eList,exaLong_t,offsetof(elm_rcb,id));

  exaBarrier(h);

  data=exaArrayGetPointer(eList);
  nel=exaArrayGetSize(eList);
  for(e=0;e<nel;e++) {
    part[e]=data[e].orig;
  }

  buffer_free(&buf);

  exaArrayFree(eList);
  exaFinalize(h);

  return 0;
}
