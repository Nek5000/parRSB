#include <genmap-impl.h>

struct csrMat_{
  GenmapUInt rn,cn,*rowOffsets,*colIdx;
  GenmapScalar *v;
};
typedef struct csrMat_ * csrMat;

static void applyCsrMat(GenmapVector y,csrMat M,GenmapVector x){
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

struct mgData_{
  GenmapUInt *levelOffsets;
  GenmapScalar *x,*b;
};
typedef struct mgData_ * mgData;

struct mgLevel_{
  mgData data;
  int nSmooth;
  GenmapUInt nLocal;
  GenmapScalar sigma;
  GenmapComm comm;
};
typedef struct mgLevel_ * mgLevel;
