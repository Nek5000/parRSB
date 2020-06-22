#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

#include <exa.h>

typedef struct {
  GenmapULong r,c; /* 1-index */
  GenmapInt owner; /* 0-index */
  GenmapInt dest; /* 0-index */
} entry;

void setDestination(struct array *entries,size_t inOffset,size_t outOffset,
  GenmapHandle h,GenmapComm c)
{
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

void setupParallelCSR(GenmapHandle h,GenmapComm c,parMat *mat)
{
  GenmapInt lelt=GenmapGetNLocalElements(h);
  GenmapInt nv  =GenmapGetNVertices(h);

  GenmapScan(h,c);
  GenmapLong startId=GenmapGetLocalStartIndex(h)+1;

  GenmapLong *eIds; GenmapInt *offsets;
  GenmapFindNeighbors(h,c,&eIds,&offsets);

  struct array entries=null_array;
  array_init(entry,&entries,offsets[lelt]);

  entry *ptr=entries.ptr;
  GenmapInt i,j,n=0;
  for(i=0;i<lelt;i++)
    for(j=0;j<offsets[i];j++)
      ptr[n].r=i+startId,ptr[n].c=((eIds[n]<0)?-eIds[n]:eIds[n]),n++;

  setDestination(&entries,offsetof(entry,r),offsetof(entry,dest ),h,c);
  setDestination(&entries,offsetof(entry,c),offsetof(entry,owner),h,c);

  GenmapFree(eIds); GenmapFree(offsets);
  array_free(&entries);
}

static void applyMat(GenmapVector x,parMat M,GenmapVector y){
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
