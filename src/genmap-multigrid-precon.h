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
  buffer buf;
};

void parMatSetup(GenmapHandle h,GenmapComm c,parMat *M);
void parMatApply(GenmapVector x,parMat M,GenmapVector y);
void parMatPrint(parMat M);
int  parMatFree(parMat M);

struct mgLevel_{
  int nSmooth;
  GenmapScalar sigma;
  struct comm c;
  parMat M;
};

struct mgData_{
  GenmapInt nLevels;
  mgLevel *levels;
  GenmapUInt *levelOffsets;
  GenmapScalar *x,*b;
};

void mgSetup(GenmapComm c,parMat M,mgData *d);
void mgLevelSetup(mgLevel l0,mgLevel *l1_);
void mgFree(mgData d);

int log2i(GenmapInt i);

#endif
