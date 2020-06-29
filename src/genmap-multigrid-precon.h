#ifndef _GENMAP_PRECON_H_
#define _GENMAP_PRECON_H_

typedef struct parMat_  *parMat;
typedef struct mgData_  *mgData;
typedef struct mgLevel_ *mgLevel;

struct parMat_{
  GenmapUInt rn;
  GenmapULong cn;

  GenmapULong rowStart;
  sint rank;

  GenmapUInt *rowOffsets;
  GenmapULong *colIdx;

  sint *owner;
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
  struct gs_data *J; // interpolation from level i to i+1
  parMat M;
};

struct mgData_{
  struct comm c;
  sint nLevels;
  mgLevel *levels;
  GenmapUInt *levelOffsets;
  GenmapScalar *x,*b,*buf;
};

void mgSetup(GenmapComm c,parMat M,mgData *d);
void mgLevelSetup(mgLevel l0,mgLevel *l1_);
void mgFree(mgData d);

int log2i(sint i);

typedef struct{
  GenmapULong r,c,rn,cn;
  sint owner; /* 0-index */
  sint p;
  GenmapScalar v;
} entry;

#define GETLNG(ptr,i,offset)\
  (*((GenmapULong*)((char*)(ptr)+offset+(i)*sizeof(entry))))

#define GETPTR(ptr,i,offset) ((char*)(ptr)+offset+i*sizeof(entry))

void setOwner(char *ptr,sint n,size_t inOffset,size_t outOffset,
  GenmapLong lelg,sint np);

#endif
