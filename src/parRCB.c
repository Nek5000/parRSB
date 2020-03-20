#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "exa-impl.h"
#include "parRSB.h"

void fparRCB_partMesh(int *part,double *vtx,int *nel,int *ndim,
  int *options,int *comm,int *err)
{
  *err = 1;

  GenmapCommExternal c;
  c = MPI_Comm_f2c(*comm);
  *err = parRCB_partMesh(part,vtx,*nel,*ndim,options,c);
}

int parRCB_partMesh(int *part,double *vtx,int nel,int ndim,
  int *options,MPI_Comm comm)
{
  exaHandle h;
  exaInit(&h,comm,"/host");

  int rank=exaRank(h),size=exaSize(h);

  /* load balance input data */
  exaLong nelg;
  exaLong nell = nel;
  exaReduce(h,&nelg,&nell,1,exaLong_t,exaAddOp);
  exaLong nstar = nelg/size;
  if(nstar == 0) nstar = 1;

  exaLong nelg_start,buf0;
  exaScan(h,&nelg_start,&nell,&buf0,1,exaLong_t,exaAddOp);
  nelg_start-=nell;

  struct array eList;
  elm_rcb *data;

  array_init(elm_rcb,&eList,nel), eList.n=nel;
  int e, n;
  for(data=eList.ptr,e=0; e<nel; ++e) {
    data[e].id=nelg_start+(e+1);
    const exaLong eg=data[e].id;

    data[e].proc=(int) ((eg-1)/nstar);
    if(eg>size*nstar) data[e].proc= (eg%size)-1; 
    for(int n=0;n<ndim;n++)
      data[e].coord[n]=vtx[e*ndim+n];
  }

  struct comm c; comm_init(&c,comm);
  struct crystal cr; crystal_init(&cr,&c);
  buffer buf; buffer_init(&buf,1024);

  sarray_transfer(elm_rcb,&eList,proc,1,&cr);
  data=eList.ptr;
  nel =eList.n;

  sarray_sort(elm_rcb,data,(unsigned int)nel,id,TYPE_LONG,&buf);

  MPI_Comm commRSB;
  MPI_Comm_split(comm, nel>0, rank, &commRSB);

  if(nel>0) {
    double time0 = comm_time();
    if(options!=NULL) exaSetDebug(h,options[0]);
    if(exaRank(h)==0)
      exaDebug(h,"\nfinished in %lfs\n",comm_time()-time0);

    fflush(stdout);
  }

  MPI_Comm_free(&commRSB);

  /* restore original input */
  sarray_transfer(elm_rcb,&eList,proc,0,&cr);
  data=eList.ptr;
  nel =eList.n;
  sarray_sort(elm_rcb,data,(unsigned int)nel,id,TYPE_LONG,&buf);
  MPI_Barrier(comm);

  for(e = 0; e < nel; e++) {
    part[e]=data[e].part;
  }

  array_free(&eList);
  buffer_free(&buf);
  crystal_free(&cr);
  comm_free(&c);
  exaFinalize(h);

  return 0;
}
