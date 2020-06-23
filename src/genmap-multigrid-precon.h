#ifndef _GENMAP_PRECON_H_
#define _GENMAP_PRECON_H_

typedef struct parMat_  *parMat;
typedef struct mgData_  *mgData;
typedef struct mgLevel_ *mgLevel;

struct parMat_{
  GenmapUInt rn;
  GenmapULong cn;

  GenmapULong rowStart;
  GenmapInt rank;

  GenmapUInt *rowOffsets;
  GenmapULong *colIdx;

  GenmapInt *owner;
  GenmapScalar *v,*x;

  struct gs_data *gsh;
  GenmapComm c;
  GenmapHandle h;

  buffer buf;
};

void parMatSetup(GenmapHandle h,GenmapComm c,parMat *M);
void parMatApply(GenmapVector x,parMat M,GenmapVector y);
void parMatPrint(parMat M);
int parMatFree(parMat M);

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
