#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

#include <exa.h>

#define ABS(i) ((i<0)?-i:i)

void setOwner(char *ptr,GenmapInt n,size_t inOffset,size_t outOffset,
  GenmapLong lelg,GenmapInt np)
{
  GenmapInt lelt=lelg/np;
  GenmapInt nrem=lelg%np;

  GenmapULong *inPtr;
  GenmapInt  *outPtr;
  GenmapInt i; GenmapLong row;
  for(i=0; i<n; i++){
    inPtr =(GenmapULong*)GETPTR(ptr,i,inOffset );
    outPtr=(GenmapInt  *)GETPTR(ptr,i,outOffset);
    row   =*inPtr-1;
    //FIXME: Assumes the 'reverse-Nek' element distribution
#if 0
    if(row<lelt*(np-nrem)) *outPtr=(GenmapInt) row/lelt;
    else *outPtr=np-nrem+(GenmapInt) (row-lelt*(np-nrem))/(lelt+1);
#else
    if(nrem==0) *outPtr=(GenmapInt) row/lelt;
    else if(row<(lelt+1)*nrem) *outPtr=(GenmapInt) row/(lelt+1);
    else *outPtr=nrem+(GenmapInt) (row-(lelt+1)*nrem)/lelt;
#endif
  }
}

void parMatSetup(GenmapHandle h,GenmapComm c,parMat *M_)
{
  GenmapInt lelt=GenmapGetNLocalElements(h);
  GenmapInt nv  =GenmapGetNVertices(h);

  GenmapScan(h,c);
  GenmapLong s=GenmapGetLocalStartIndex(h)+1;

  GenmapLong *eIds; GenmapInt *offsets;
  GenmapFindNeighbors(h,c,&eIds,&offsets);

  struct array entries=null_array;
  array_init(entry,&entries,offsets[lelt]);
  entries.n=offsets[lelt];

  entry *ptr=entries.ptr;
  GenmapInt i,j,n=0;
  for(i=0;i<lelt;i++)
    for(j=0;j<offsets[i]+1;j++){
      ptr[n].r=s+i,ptr[n].c=ABS(eIds[n]),ptr[n].v=-1.0;
      if(ptr[n].r==ptr[n].c) ptr[n].v=offsets[i];
      n++;
    }
  GenmapFree(offsets);

  GenmapScan(h,c);
  GenmapLong ng=GenmapGetNGlobalElements(h);
  GenmapInt np =GenmapCommSize(c);
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

  M->rn=n,M->cn=GenmapGetNGlobalElements(h);
  M->rank=GenmapCommRank(c);

  GenmapScan(h,c);
  M->rowStart=GenmapGetLocalStartIndex(h)+1;

  GenmapMalloc(M->rn+1,&M->rowOffsets);

  if(n==0)
    M->colIdx=NULL,M->v=NULL,M->owner=NULL;
  else
    GenmapMalloc(entries.n,&M->colIdx),GenmapMalloc(entries.n,&M->v),\
      GenmapMalloc(entries.n,&M->owner );

  ptr=entries.ptr;
  for(i=0; i<entries.n; i++)
    M->colIdx[i]=ptr[i].c,M->v[i]=ptr[i].v,M->owner[i]=ptr[i].owner;

  M->rowOffsets[0]=0,i=0; GenmapInt nn=0;
  while(i<entries.n){
    j=i+1;
    while(j<entries.n && ptr[i].r==ptr[j].r) j++;
    i=M->rowOffsets[++nn]=j;
  }
  assert(n==nn);
  assert(M->rowOffsets[n]=entries.n);

  if(entries.n>0) GenmapRealloc(entries.n,&eIds);

  for(i=0; i<M->rn; i++)
    for(j=M->rowOffsets[i]; j<M->rowOffsets[i+1]; j++)
      if(M->rowStart+i==M->colIdx[j]) eIds[j]=M->colIdx[j];
      else eIds[j]=-M->colIdx[j];

  M->gsh=gs_setup(eIds,M->rowOffsets[n],&c->gsc,0,gs_crystal_router,0);
  buffer_init(&M->buf,1024);

  GenmapFree(eIds);
  array_free(&entries);
}

void parMatApply(GenmapScalar *y,parMat M,GenmapScalar *x,
  GenmapScalar *buf)
{
  const GenmapUInt rn=M->rn,cn=M->cn;
  const GenmapUInt *offsets=M->rowOffsets;
  const GenmapScalar *v=M->v;

  GenmapLong s=M->rowStart;
  GenmapInt i,j;
  for(i=0;i<M->rn;i++)
    for(j=M->rowOffsets[i]; j<M->rowOffsets[i+1]; j++)
      if(M->owner[j]==M->rank && M->colIdx[j]==s+i) buf[j]=x[i];

#if defined(EXA_DEBUG)
  printf("buf(before): ");
  for(i=0;i<M->rowOffsets[rn];i++)
    printf("%lf ",buf[i]);
  printf("\n");
#endif

  gs(buf,genmap_gs_scalar,gs_add,0,M->gsh,&M->buf);

#if defined(EXA_DEBUG)
  printf("buf(after): ");
  for(i=0;i<M->rowOffsets[rn];i++)
    printf("%lf ",buf[i]);
  printf("\n");
#endif

  GenmapUInt je;
  for(i=0;i<rn;i++){
    for(y[i]=0.0,j=offsets[i],je=offsets[i+1]; j<je; j++)
      y[i]+=(*v++)*(*buf++);
  }
}

void parMatPrint(parMat M){
  const GenmapInt rn=M->rn;
  const GenmapUInt *offsets=M->rowOffsets;
  const GenmapScalar *v=M->v;
  const GenmapULong *col=M->colIdx;

  printf("rn=%d\n",M->rn);
  GenmapInt i,j;
  for(i=0; i<rn; i++){
    printf("row %d: ",i+1);
    for(j=offsets[i]; j<offsets[i+1]; j++)
      printf("(%lld,%lf)",col[j],v[j]);
    printf("\n");
  }
}

int parMatFree(parMat M){
  if(M->colIdx) GenmapFree(M->colIdx);
  if(M->v) GenmapFree(M->v);
  if(M->rowOffsets) GenmapFree(M->rowOffsets);
  if(M->owner) GenmapFree(M->owner);
  if(M->gsh) gs_free(M->gsh);
  if(M->buf.ptr) buffer_free(&M->buf);
  GenmapFree(M);

  return 0;
}
