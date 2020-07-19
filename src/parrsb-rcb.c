#include "exasort-impl.h"
#include "parRSB.h"

void get_axis_len(double *length,struct array *a,struct comm *c,int ndim)
{
  double min[MAXDIM],max[MAXDIM];
  sint i;
  for(i=0;i<ndim;i++) min[i]=DBL_MAX,max[i]=-DBL_MAX;

  sint nel=a->n;
  elm_rcb* elems=a->ptr;
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

  double buf[MAXDIM];
  comm_allreduce(c,gs_double,gs_min,min,MAXDIM,buf);
  comm_allreduce(c,gs_double,gs_max,max,MAXDIM,buf);

  for(i=0;i<ndim;i++)
    length[i]=max[i]-min[i];
}

int parRCB(struct comm *ci,struct array *a,int ndim){
  struct comm c; comm_dup(&c,ci);

  uint offsets[3]={offsetof(elm_rcb,coord[0]),
    offsetof(elm_rcb,coord[1]),offsetof(elm_rcb,coord[2])};

  double length[MAXDIM];

  sint rank=c.id;
  sint size=c.np;

  while(size>1){
    get_axis_len(length,a,&c,ndim);

    int axis1=0,d;
    for(d=1;d<ndim;d++)
      if(length[d]>length[axis1]) axis1=d;
    int axis2=(axis1+1)%2;
    for(d=0;d<ndim;d++)
      if(length[d]>length[axis2] && d!=axis1) axis2=d;

    uint off=offsets[axis1];
#if 0 //FIXME
    parallel_sort(elm_rcb,a,off,gs_double,&c);
#else
    sort_data_private sd;

    sd.nfields=1;
    sd.unit_size=sizeof(elm_rcb);
    sd.align    =ALIGNOF(elm_rcb);
    sd.t[0]     =gs_double;
    sd.offset[0]=off;

    sd.a=a;

    sd.balance=1;
    sd.algo=exaSortAlgoBinSort;

    buffer_init(&sd.buf,1024);
    sort_private(&sd,&c);
    buffer_free(&sd.buf);
#endif

    int p=(size+1)/2;
    int bin=(rank>=p);

    comm_ext old=c.c;
#ifdef MPI
    MPI_Comm new; MPI_Comm_split(old,bin,rank,&new);
    comm_free(&c); comm_init(&c,new);
    MPI_Comm_free(&new);
#endif
    rank=c.id;
    size=c.np;
  }

  comm_free(&c);
  return 0;
}

void get_axis_len_local(double *length,elm_rcb *elems,uint nel,int ndim)
{
  double min[MAXDIM],max[MAXDIM];
  sint i;
  for(i=0;i<ndim;i++) min[i]=DBL_MAX,max[i]=-DBL_MAX;

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

  for(i=0;i<ndim;i++)
    length[i]=max[i]-min[i];
}

void rcb_local(struct array *a,uint start,uint end,int ndim,buffer *buf){
  assert(end>=start);

  uint size=end-start;
  if(size<=2) return;

  double length[3]; elm_rcb *st=a->ptr; st+=start;
  get_axis_len_local(length,st,size,ndim);

  int axis=0;
  if(fabs(length[axis])<fabs(length[1])) axis=1;
  if(ndim==3)
    if(fabs(length[axis])<fabs(length[2])) axis=2;

  switch(axis){
    case 0:
      sarray_sort(elm_rcb,st,size,coord[0],3,buf);
      break;
    case 1:
      sarray_sort(elm_rcb,st,size,coord[1],3,buf);
      break;
    case 2:
      sarray_sort(elm_rcb,st,size,coord[2],3,buf);
      break;
    default:
      break;
  }

  uint mid=(start+end)/2;
  rcb_local(a,start,mid,ndim,buf);
  rcb_local(a,mid  ,end,ndim,buf);
}
