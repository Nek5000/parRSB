#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

#include <exa.h>

#define ABS(i) ((i<0)?-i:i)

void setOwner(char *ptr,sint n,size_t inOffset,size_t outOffset,
  slong lelg,sint np)
{
  sint lelt=lelg/np;
  sint nrem=lelg%np;

  ulong *inPtr;
  sint  *outPtr;
  sint i; slong row;
  for(i=0; i<n; i++){
    inPtr =(ulong*)GETPTR(ptr,i,inOffset );
    outPtr=(sint  *)GETPTR(ptr,i,outOffset);
    row   =*inPtr-1;
    //FIXME: Assumes the 'reverse-Nek' element distribution
#if 0
    if(row<lelt*(np-nrem)) *outPtr=(sint) row/lelt;
    else *outPtr=np-nrem+(sint) (row-lelt*(np-nrem))/(lelt+1);
#else
    if(nrem==0) *outPtr=(sint) row/lelt;
    else if(row<(lelt+1)*nrem) *outPtr=(sint) row/(lelt+1);
    else *outPtr=nrem+(sint) (row-(lelt+1)*nrem)/lelt;
#endif
  }
}

void parMatSetup(GenmapHandle h,GenmapComm c,parMat *M_)
{
  sint lelt=GenmapGetNLocalElements(h);
  sint nv  =GenmapGetNVertices(h);

  GenmapScan(h,c);
  slong s=GenmapGetLocalStartIndex(h)+1;

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

  GenmapScan(h,c);
  slong ng=GenmapGetNGlobalElements(h);
  sint np =GenmapCommSize(c);
  setOwner(entries.ptr,entries.n,offsetof(entry,c),
    offsetof(entry,owner),ng,np);

#if defined(EXA_DEBUG)
  for(i=0; i<entries.n; i++)
    printf("elem: (%lld,%lld,%d)\n",ptr[i].r,ptr[i].c,ptr[i].owner);
#endif

  buffer buf; buffer_init(&buf,1024);
  ptr=entries.ptr;
  sarray_sort_2(entry,ptr,entries.n,r,1,c,1,&buf);
  buffer_free(&buf);

  GenmapMalloc(1,M_); parMat M=*M_;

  i=0,n=0;
  while(i<entries.n){
    j=i+1;
    while(j<entries.n && ptr[i].r==ptr[j].r) j++;
    i=j,n++;
  }

  M->rn=n;
  GenmapScan(h,c);
  M->row_start=GenmapGetLocalStartIndex(h)+1;

  GenmapMalloc(M->rn+1,&M->row_off);

  if(n==0)
    M->col=NULL,M->v=NULL;
  else
    GenmapMalloc(entries.n,&M->col),GenmapMalloc(entries.n,&M->v);

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
      if(M->row_start+i==M->col[j]) eIds[j]=M->col[j];
      else eIds[j]=-M->col[j];

  M->gsh=gs_setup(eIds,M->row_off[n],&c->gsc,0,gs_crystal_router,0);
  buffer_init(&M->buf,1024);

  GenmapFree(eIds);
  array_free(&entries);
}

void parMatApply(GenmapScalar *y,parMat M,GenmapScalar *x,
  GenmapScalar *buf)
{
  const uint rn=M->rn;
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

void parMatPrint(parMat M,struct comm *c){
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

int parMatFree(parMat M){
  if(M->col) GenmapFree(M->col);
  if(M->v) GenmapFree(M->v);
  if(M->row_off) GenmapFree(M->row_off);
  if(M->gsh) gs_free(M->gsh);
  if(M->buf.ptr) buffer_free(&M->buf);
  GenmapFree(M);

  return 0;
}
