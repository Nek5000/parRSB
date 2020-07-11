#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

#include <exa.h>

#define ABS(i) ((i<0)?-i:i)

void csr_mat_setup(GenmapHandle h,GenmapComm c,csr_mat *M_)
{
  sint lelt=GenmapGetNLocalElements(h);
  sint nv  =GenmapGetNVertices(h);

  GenmapScan(h,c); slong s=GenmapGetLocalStartIndex(h)+1;

  slong *eIds; sint *offsets;
  GenmapFindNeighbors(h,c,&eIds,&offsets);

  struct array entries=null_array;
  array_init(entry,&entries,offsets[lelt]);
  entries.n=offsets[lelt];

  entry *ptr=entries.ptr;
  sint i,j,n=0;
  for(i=0;i<lelt;i++)
    for(j=0;j<offsets[i]+1;j++){
      ptr[n].r=s+i,ptr[n].c=ABS(eIds[n]),ptr[n].v=-1.0;
      if(ptr[n].r==ptr[n].c) ptr[n].v=offsets[i];
      n++;
    }
  GenmapFree(offsets);

  GenmapScan(h,c); slong ng=GenmapGetNGlobalElements(h);
  sint np=GenmapCommSize(c);

  buffer buf; buffer_init(&buf,1024);
  ptr=entries.ptr;
  sarray_sort_2(entry,ptr,entries.n,r,1,c,1,&buf);
  buffer_free(&buf);

  GenmapMalloc(1,M_); csr_mat M=*M_;

  i=0,n=0;
  while(i<entries.n){
    j=i+1;
    while(j<entries.n && ptr[i].r==ptr[j].r) j++;
    i=j,n++;
  }

  M->rn=n;
  GenmapScan(h,c); M->row_start=GenmapGetLocalStartIndex(h)+1;

  GenmapMalloc(M->rn+1,&M->row_off);

  if(n==0)
    M->col=NULL,M->v=NULL,M->diag=NULL;
  else{
    GenmapMalloc(entries.n,&M->col),GenmapMalloc(entries.n,&M->v);
    GenmapMalloc(M->rn,&M->diag);
  }

  ptr=entries.ptr;
  for(i=0; i<entries.n; i++)
    M->col[i]=ptr[i].c,M->v[i]=ptr[i].v;

  M->row_off[0]=0,i=0; sint nn=0;
  while(i<entries.n){
    j=i+1;
    while(j<entries.n && ptr[i].r==ptr[j].r) j++;
    i=M->row_off[++nn]=j;
  }
  assert(n==nn);
  assert(M->row_off[n]=entries.n);

  if(entries.n>0) GenmapRealloc(entries.n,&eIds);

  for(i=0; i<M->rn; i++)
    for(j=M->row_off[i]; j<M->row_off[i+1]; j++)
      if(M->row_start+i==M->col[j])
        eIds[j]=M->col[j],M->diag[i]=M->v[j];
      else eIds[j]=-M->col[j];

  M->gsh=gs_setup(eIds,M->row_off[n],&c->gsc,0,gs_crystal_router,0);
  buffer_init(&M->buf,1024);

  GenmapFree(eIds);
  array_free(&entries);
}

struct gs_data *get_csr_top(csr_mat M,struct comm *c){
  const uint rn=M->rn;
  const uint n =M->row_off[rn];

  slong *ids;
  if(n>0) GenmapMalloc(n,&ids);

  uint i,j;
  for(i=0; i<rn; i++)
    for(j=M->row_off[i]; j<M->row_off[i+1]; j++)
      if(M->row_start+i==M->col[j]) ids[j]=M->col[j];
      else ids[j]=-M->col[j];

  struct gs_data *gsh=gs_setup(ids,n,c,0,gs_auto,0);

  if(n>0) GenmapFree(ids);

  return gsh;
}

void csr_mat_apply(GenmapScalar *y,csr_mat M,GenmapScalar *x,
  GenmapScalar *buf)
{
  const uint rn=M->rn;
  if(rn==0) return;

  const uint *offsets=M->row_off;
  const GenmapScalar *v=M->v;

  ulong s=M->row_start;
  sint i,j;
  for(i=0;i<M->rn;i++)
    for(j=M->row_off[i]; j<M->row_off[i+1]; j++)
      if(M->col[j]==s+i) buf[j]=x[i];

#if defined(EXA_DEBUG)
  printf("buf(before): ");
  for(i=0;i<M->row_off[rn];i++)
    printf("%lf ",buf[i]);
  printf("\n");
#endif

  gs(buf,genmap_gs_scalar,gs_add,0,M->gsh,&M->buf);

#if defined(EXA_DEBUG)
  printf("buf(after): ");
  for(i=0;i<M->row_off[rn];i++)
    printf("%lf ",buf[i]);
  printf("\n");
#endif

  uint je;
  for(i=0;i<rn;i++){
    for(y[i]=0.0,j=offsets[i],je=offsets[i+1]; j<je; j++)
      y[i]+=(*v++)*(*buf++);
  }
}

void csr_mat_print(csr_mat M,struct comm *c){
  const sint rn=M->rn;
  const uint *offsets=M->row_off;
  const GenmapScalar *v=M->v;
  const ulong *col=M->col;

  uint i,j,k;

  for(k=0; k<c->np; k++){
    comm_barrier(c);
    if(c->id==k)
      for(i=0; i<rn; i++){
        printf("%02u: r%02lld: ",c->id,M->row_start+i);
        for(j=offsets[i]; j<offsets[i+1]; j++)
          printf("(%2lld,%.2lf)",col[j],v[j]);
        printf("\n");
      }
  }
}

int csr_mat_free(csr_mat M){
  if(M->col) GenmapFree(M->col);
  if(M->v) GenmapFree(M->v);
  if(M->diag) GenmapFree(M->diag);
  if(M->row_off) GenmapFree(M->row_off);
  if(M->gsh) gs_free(M->gsh);
  if(M->buf.ptr) buffer_free(&M->buf);
  GenmapFree(M);

  return 0;
}
