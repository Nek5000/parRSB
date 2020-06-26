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
  GenmapScalar *v;

  struct gs_data *gsh;
  buffer buf;
};

void parMatSetup(GenmapHandle h,GenmapComm c,parMat *M);
void parMatApply(GenmapScalar *y,parMat M,GenmapScalar *x,
  GenmapScalar *buf);
void parMatPrint(parMat M);
int  parMatFree(parMat M);

struct mgLevel_{
  mgData data;
  int nSmooth;
  GenmapScalar sigma;
  parMat M;
};

struct mgData_{
  struct comm c;
  GenmapInt nLevels;
  mgLevel *levels;
  GenmapUInt *levelOffsets;
  GenmapScalar *x,*b;
};

void mgSetup(GenmapComm c,parMat M,mgData *d);
void mgLevelSetup(mgLevel l0,mgLevel *l1_);
void mgFree(mgData d);

int log2i(GenmapInt i);

void setOwner(char *ptr,GenmapInt n,size_t inOffset,size_t outOffset,
  GenmapLong lelg,GenmapInt np);

#endif
