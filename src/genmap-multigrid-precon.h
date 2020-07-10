#ifndef _GENMAP_PRECON_H_
#define _GENMAP_PRECON_H_

typedef struct parMat_  *parMat;
typedef struct mgData_  *mgData;
typedef struct mgLevel_ *mgLevel;

struct parMat_{
  uint rn;
  ulong row_start;

  uint *row_off;
  ulong *col;
  GenmapScalar *v,*diag;

  struct gs_data *gsh;
  buffer buf;
};

// for the coarse level
void parMatSetup(GenmapHandle h,GenmapComm c,parMat *M);
void parMatApply(GenmapScalar *y,parMat M,GenmapScalar *x,
  GenmapScalar *buf);
void parMatPrint(parMat M,struct comm *c);
int  parMatFree(parMat M);

struct mgLevel_{
  mgData data;
  int nsmooth;
  GenmapScalar sigma;
  struct gs_data *J; // interpolation from level i to i+1
  parMat M;
};

struct mgData_{
  struct comm c;
  sint nlevels;
  mgLevel *levels;
  uint *level_off;
  GenmapScalar *y,*x,*b,*buf;
};

void mgSetup(GenmapComm c,parMat M,mgData *d);
void mgLevelSetup(mgData data,uint level);
void mgFree(mgData d);

int log2i(sint i);

typedef struct{
  ulong r,c,rn,cn;
  uint p;
  GenmapScalar v;
} entry;

#define GETLNG(ptr,i,offset)\
  (*((ulong*)((char*)(ptr)+offset+(i)*sizeof(entry))))

#define GETPTR(ptr,i,offset) ((char*)(ptr)+offset+i*sizeof(entry))

void setOwner(char *ptr,sint n,size_t inOffset,size_t outOffset,
  slong lelg,sint np);

#endif
