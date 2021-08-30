#ifndef _GENMAP_PRECON_H_
#define _GENMAP_PRECON_H_

#include <genmap.h>

typedef struct mgData_ *mgData;
typedef struct mgLevel_ *mgLevel;

struct csr_mat {
  uint rn;
  ulong row_start;

  uint *roff;
  ulong *col;
  GenmapScalar *v, *diag, *buf;

  struct gs_data *gsh;
};

// for the coarse level
void csr_mat_setup(struct csr_mat *M, struct array *entries, struct comm *c,
                   buffer *buf);
void csr_mat_apply(GenmapScalar *y, struct csr_mat *M, GenmapScalar *x,
                   buffer *buf);
void csr_mat_print(struct csr_mat *M, struct comm *c);
int csr_mat_free(struct csr_mat *M);

struct mgLevel_ {
  mgData data;
  int nsmooth;
  GenmapScalar sigma;
  struct gs_data *J; // interpolation from level i to i+1
  struct gs_data *Q; // global to local conversion of a vector
  struct csr_mat *M;
};

struct mgData_ {
  struct comm c;
  struct gs_data *top;
  buffer bfr;
  int nlevels;
  mgLevel *levels;
  uint *level_off;
  GenmapScalar *y, *x, *b, *u, *rhs, *buf;
};

void mgSetup(struct comm *c, struct csr_mat *M, mgData *d);
void mgLevelSetup(mgData data, uint level);
void mgFree(mgData d);

typedef struct {
  ulong r, c;
  uint proc;
} csr_entry;

typedef struct {
  ulong r, c, rn, cn;
  uint p;
  GenmapScalar v;
} entry;

#define GETLNG(p, i, off)                                                      \
  (*((ulong *)((char *)(p) + (off) + (i) * sizeof(entry))))
#define GETPTR(p, i, off) ((char *)(p) + (off) + (i) * sizeof(entry))

void setOwner(char *ptr, sint n, size_t inOffset, size_t outOffset, slong lelg,
              sint np);

void mg_vcycle(GenmapScalar *u, GenmapScalar *rhs, mgData d);

#endif
