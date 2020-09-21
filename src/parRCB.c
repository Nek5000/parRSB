#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include <sort.h>
#include <genmap-impl.h>
#include <parRSB.h>

void fparRCB_partMesh(int *part,int *seq,double *vtx,int *nel,int *nv,
  int *options,int *comm,int *err)
{
  *err = 1;

  comm_ext c; c = MPI_Comm_f2c(*comm);
  *err=parRCB_partMesh(part,seq,vtx,*nel,*nv,options,c);
}

// vtx = [nel,nv,ndim]
int parRCB_partMesh(int *part,int *seq,double *vtx,int nel,int nv,
  int *options,MPI_Comm comm)
{
  struct comm c; comm_init(&c,comm);
  int rank=c.id,size=c.np;

  slong out[2][1],buf[2][1];
  slong nell=nel;
  comm_scan(out,&c,gs_long,gs_add,&nell,1,&buf);
  slong nelg_start=out[0][0];
  slong nelg      =out[1][0];

  struct array a; array_init(elm_rcb,&a,nel);
  elm_rcb *data=a.ptr;

  int ndim=(nv==8)?3:2;

  int e,n,v;
  for(e=0;e<nel;++e){
    data[e].id=nelg_start+(e+1);
    data[e].orig=rank;
    data[e].coord[0]=data[e].coord[1]=data[e].coord[2]=0.0;
    for(v=0;v<nv;v++){
      for(n=0;n<ndim;n++)
        data[e].coord[n]+=vtx[e*ndim*nv+v*ndim+n];
    }
    for(n=0;n<ndim;n++)
      data[e].coord[n]/=nv;
  }
  a.n=nel;

  //TODO: load balance

  buffer bfr; buffer_init(&bfr,1024);
  metric_init();

  struct comm rcb;
  comm_ext old=c.c;
#ifdef MPI
  MPI_Comm new; MPI_Comm_split(old,nel>0,rank,&new);
  comm_init(&rcb,new); MPI_Comm_free(&new);
#else
  comm_init(&rcb,1);
#endif

  if(nel>0){
    comm_barrier(&rcb);
    double time=comm_time();

    if(c.id==0)
      printf("running RCB ...");
    fflush(stdout);

    parRCB(&rcb,&a,ndim);

    // do local rcb
    uint s1=0,e1=a.n;
    rcb_local(&a,s1,e1,ndim,&bfr);

    elm_rcb *ptr=a.ptr;
    int i;
    for(i=0; i<a.n; i++)
      ptr[i].seq=i;

    comm_barrier(&rcb);
    time=comm_time()-time;

    if(c.id==0)
      printf(" finished in %g s\n",time);
    fflush(stdout);
  }

  /* restore original input */
  struct crystal cr; crystal_init(&cr,&c);
  sarray_transfer(elm_rcb,&a,orig,1,&cr);
  crystal_free(&cr);
  assert(a.n==nel);
  sarray_sort(elm_rcb,a.ptr,a.n,id,1,&bfr);


  data=a.ptr;
  for(e=0;e<nel;e++){
    part[e]=data[e].orig;
    seq [e]=data[e].seq ;
  }

  metric_print(&c);

  comm_free(&rcb);
  metric_finalize();
  buffer_free(&bfr);
  comm_free(&c);

  return 0;
}
