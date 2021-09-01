#ifndef _GENMAP_PRECON_H_
#define _GENMAP_PRECON_H_

#include <genmap-impl.h>

typedef struct mgData_ *mgData;
typedef struct mgLevel_ *mgLevel;

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

#define GETLNG(p, i, off)                                                      \
  (*((ulong *)((char *)(p) + (off) + (i) * sizeof(entry))))
#define GETPTR(p, i, off) ((char *)(p) + (off) + (i) * sizeof(entry))

void setOwner(char *ptr, sint n, size_t inOffset, size_t outOffset, slong lelg,
              sint np);

void mg_vcycle(GenmapScalar *u, GenmapScalar *rhs, mgData d);

/* Iterative methods */
int flex_cg(genmap_vector x, struct gs_laplacian *gl, mgData d,
            genmap_vector ri, int maxIter, struct comm *gsc, buffer *buf);

int project(genmap_vector x, struct gs_laplacian *gl, mgData d,
            genmap_vector ri, int max_iter, struct comm *gsc, buffer *buf);

#endif
