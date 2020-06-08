#include <genmap-impl.h>

struct csrMat_{
  GenmapUInt rn,cn,*rowOffsets,*colIdx;
  double *v;
};
typedef struct csrMat_ * csrMat;

static void applyM(GenmapVector y,csrMat M,GenmapVector x){
  GenmapUInt i; const GenmapUInt rn=M->rn,cn=M->cn;
  const GenmapUInt *offsets=M->rowOffsets,*col=M->colIdx;
}
