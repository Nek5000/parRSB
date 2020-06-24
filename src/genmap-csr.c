#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

#include <exa.h>

#define ABS(i) ((i<0)?-i:i)

typedef struct {
  GenmapULong r,c; /* 1-index */
  GenmapInt owner; /* 0-index */
  GenmapInt dest; /* 0-index */
  GenmapScalar v;
} entry;

#define GETPTR(ptr,i,offset) ((char*)(ptr)+offset+i*sizeof(entry))

void setDestination(struct array *entries,size_t inOffset,
  size_t outOffset,GenmapHandle h,GenmapComm c){
  GenmapScan(h,c);
  GenmapLong lelg=GenmapGetNGlobalElements(h);
  GenmapInt np=GenmapCommSize(c);

  GenmapInt lelt=lelg/np;
  GenmapInt nrem=lelg%np;

  GenmapULong *inPtr;
  GenmapInt  *outPtr;
  GenmapInt i; GenmapLong row;
  for(i=0; i<entries->n; i++){
    inPtr =(GenmapULong*)GETPTR(entries->ptr,i,inOffset );
    outPtr=(GenmapInt  *)GETPTR(entries->ptr,i,outOffset);
    row   =*inPtr-1;
    if(row<lelt*(np-nrem)) *outPtr=(GenmapInt) row/lelt;
    else *outPtr=np-nrem+(GenmapInt) (row-lelt*(np-nrem))/(lelt+1);
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

  setDestination(&entries,offsetof(entry,r),offsetof(entry,dest ),h,c);
  setDestination(&entries,offsetof(entry,c),offsetof(entry,owner),h,c);

#if defined(EXA_DEBUG)
  for(i=0; i<entries.n; i++)
    printf("elem: (%lld,%lld,%d,%d)\n",ptr[i].r,ptr[i].c,ptr[i].dest,
        ptr[i].owner);
#endif

  struct crystal cr; crystal_init(&cr,&c->gsc);
  sarray_transfer(entry,&entries,dest,1,&cr);
  crystal_free(&cr);

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

  if(n==0) M->colIdx=NULL,M->v=NULL,M->x=NULL,M->owner=NULL;
  else{
    GenmapMalloc(entries.n,&M->colIdx);
    GenmapMalloc(entries.n,&M->v     );
    GenmapMalloc(entries.n,&M->x     );
    GenmapMalloc(entries.n,&M->owner );
  }

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

  for(i=0; i<entries.n; i++)
    if(ptr[i].r==ptr[i].c && ptr[i].owner==M->rank) eIds[i]=ptr[i].c;
    else eIds[i]=-ptr[i].c;

  M->gsh=gs_setup(eIds,M->rowOffsets[n],&c->gsc,0,gs_crystal_router,0);
  buffer_init(&M->buf,1024);

  GenmapFree(eIds);
  array_free(&entries);
}

void parMatApply(GenmapVector x,parMat M,GenmapVector y){
  const GenmapUInt rn=M->rn,cn=M->cn;
  const GenmapUInt *offsets=M->rowOffsets;
  const GenmapScalar *v=M->v;
  GenmapScalar *xx=M->x;

  GenmapScalar *xd=x->data;
  GenmapScalar *yd=y->data;

  GenmapLong s=M->rowStart;
  GenmapInt i,j;
  for(i=0;i<M->rn;i++)
    for(j=M->rowOffsets[i]; j<M->rowOffsets[i+1]; j++)
      if(M->owner[j]==M->rank && M->colIdx[j]==s+i) xx[j]=xd[i];

#if defined(EXA_DEBUG)
  printf("xx(before): ");
  for(i=0;i<M->rowOffsets[rn];i++)
    printf("%lf ",xx[i]);
  printf("\n");
#endif

  gs(xx,genmap_gs_scalar,gs_add,0,M->gsh,&M->buf);

#if defined(EXA_DEBUG)
  printf("xx(after): ");
  for(i=0;i<M->rowOffsets[rn];i++)
    printf("%lf ",xx[i]);
  printf("\n");
#endif

  GenmapUInt je;
  for(i=0;i<rn;i++){
    for(yd[i]=0.0,j=offsets[i],je=offsets[i+1]; j<je; j++)
      yd[i]+=(*v++)*(*xx++);
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
  if(M->x) GenmapFree(M->x);
  if(M->rowOffsets) GenmapFree(M->rowOffsets);
  if(M->owner) GenmapFree(M->owner);
  if(M->gsh) gs_free(M->gsh);
  if(M->buf.ptr) buffer_free(&M->buf);
  GenmapFree(M);

  return 0;
}
