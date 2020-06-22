#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

#include <exa.h>

struct matEntry_{
  GenmapULong r,c;
  GenmapInt owner;
  GenmapInt ownerRank;
};
typedef struct matEntry_* matEntry;

void setupParallelCSR(exaArray csrEntries,parMat mat){
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
