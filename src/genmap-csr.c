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

void setDestination(struct array *entries,size_t inOffset,
  size_t outOffset,GenmapHandle h,GenmapComm c){
  GenmapScan(h,c);
  GenmapLong lelg=GenmapGetNGlobalElements(h);
  GenmapInt np=GenmapCommSize(c);

  GenmapInt lelt=lelg/np;
  GenmapInt nrem=lelg%np;

  char *inPtr=entries->ptr+inOffset,*outPtr=entries->ptr+outOffset;
  GenmapInt i;
  for(i=0; i<lelt; i++){
    GenmapLong row=*inPtr-1;
    if(row<lelt*(np-nrem))
      *outPtr=row/lelt;
    else
      *outPtr=np-nrem+(row-lelt*(np-nrem))/(lelt+1);
    inPtr+=inOffset,outPtr+=outOffset;
  }
}

void setupParallelCSR(GenmapHandle h,GenmapComm c,parMat *M_)
{
  GenmapInt lelt=GenmapGetNLocalElements(h);
  GenmapInt nv  =GenmapGetNVertices(h);

  GenmapScan(h,c);
  GenmapLong s=GenmapGetLocalStartIndex(h)+1;

  GenmapLong *eIds; GenmapInt *offsets;
  GenmapFindNeighbors(h,c,&eIds,&offsets);

  struct array entries=null_array;
  array_init(entry,&entries,offsets[lelt]);

  entry *ptr=entries.ptr;
  GenmapInt i,j,n=0;
  for(i=0;i<lelt;i++)
    for(j=0;j<offsets[i];j++){
      ptr[n].r=s+i,ptr[n].c=ABS(eIds[n]),ptr[n].v=-1.0;
      if(ptr[n].r==ptr[n].c) ptr[n].v=offsets[i];
      n++;
    }
  GenmapFree(offsets);

  setDestination(&entries,offsetof(entry,r),offsetof(entry,dest ),h,c);
  setDestination(&entries,offsetof(entry,c),offsetof(entry,owner),h,c);

  struct crystal cr; crystal_init(&cr,&c->gsc);
  sarray_transfer(entry,&entries,dest,1,&cr);
  crystal_free(&cr);

  buffer buf; buffer_init(&buf,1024);
  ptr=entries.ptr;
  sarray_sort_2(entry,ptr,entries.n,r,1,c,1,&buf);
  buffer_free(&buf);

  i=n=0;
  while(i<entries.n){
    j=i+1;
    while(j<entries.n && ptr[i].r==ptr[j].r) j++;
    i=j,n++;
  }

  GenmapMalloc(1,M_); parMat M=*M_;
  M->rn=n;
  M->cn=GenmapGetNGlobalElements(h);
  GenmapMalloc(M->rn+1  ,&M->rowOffsets);
  GenmapMalloc(entries.n,&M->colIdx    );
  GenmapMalloc(entries.n,&M->v         );

  for(i=0; i<entries.n; i++)
    M->colIdx[i]=ptr[i].c,M->v[i]=ptr[i].v;

  M->rowOffsets[0]=0,j=1;
  for(i=1; i<entries.n; i++){
    if(ptr[i-1].r!=ptr[i].r)
      M->rowOffsets[j++]=i;
  }
  assert(j==n);
  M->rowOffsets[n]=entries.n;

  for(i=0; i<entries.n; i++)
    if(ptr[i].owner==GenmapCommRank(c)) eIds[i]=ptr[i].c;
    else eIds[i]=-ptr[i].c;

  M->gsh=gs_setup(eIds,M->rowOffsets[n],&c->gsc,0,gs_crystal_router,0);

  GenmapFree(eIds);
  array_free(&entries);
}

void applyMat(GenmapVector x,parMat M,GenmapVector y){
  const GenmapUInt rn=M->rn,cn=M->cn;
  const GenmapUInt *offsets=M->rowOffsets,*col=M->colIdx;
  const GenmapScalar *v=M->v;

  GenmapScalar *xd=x->data;
  GenmapScalar *yd=y->data;

  GenmapUInt i,j,je;
  for(i=0;i<rn;i++){
    for(yd[i]=0.0,j=offsets[i],je=offsets[i+1]; j<je; j++)
      yd[i]+=(*v++)*xd[*col++];
  }
}

int freeMat(parMat M){
  GenmapFree(M->colIdx    );
  GenmapFree(M->v         );
  GenmapFree(M->rowOffsets);
  gs_free(M->gsh);
  GenmapFree(M);
}
