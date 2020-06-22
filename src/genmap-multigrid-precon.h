#ifndef _GENMAP_PRECON_H_
#define _GENMAP_PRECON_H_

typedef struct parMat_  *parMat;
typedef struct mgData_  *mgData;
typedef struct mgLevel_ *mgLevel;

struct parMat_{
  GenmapUInt rn,cn,*rowOffsets,*colIdx;
  GenmapScalar *v;
  struct gs_data *gsh;
};
static void applyMat(GenmapVector x,parMat M,GenmapVector y);

struct mgData_{
  int nLevels;
  mgLevel *levels;
  GenmapUInt *levelOffsets;
  GenmapScalar *x,*b;
};

struct mgLevel_{
  mgData data;
  int nSmooth;
  GenmapUInt nLocal;
  GenmapScalar sigma;
  GenmapComm comm;
  parMat M;
};

#endif
