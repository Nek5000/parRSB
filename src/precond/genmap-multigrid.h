#ifndef _GENMAP_PRECON_H_
#define _GENMAP_PRECON_H_

#include <genmap-impl.h>

// for the coarse level
void csr_mat_setup(struct csr_mat *M, struct array *entries, struct comm *c,
                   buffer *buf);
void csr_mat_apply(GenmapScalar *y, struct csr_mat *M, GenmapScalar *x,
                   buffer *buf);
void csr_mat_print(struct csr_mat *M, struct comm *c);
int csr_mat_free(struct csr_mat *M);

struct mg_lvl {
  int nsmooth;
  GenmapScalar sigma;
  struct gs_data *J; // interpolation from level i to i+1
  struct gs_data *Q; // global to local conversion of a vector
  struct csr_mat *M;
};
typedef struct mg_lvl *mgLevel;

struct mg_data {
  struct comm c;
  struct gs_data *top;
  buffer bfr;
  int nlevels;
  mgLevel *levels;
  uint *level_off;
  GenmapScalar *y, *x, *b, *u, *rhs, *buf;
};

void mg_setup(struct mg_data *d, int factor, struct comm *c, struct csr_mat *M);
void mg_vcycle(GenmapScalar *u, GenmapScalar *rhs, struct mg_data *d);
void mg_free(struct mg_data *d);

#endif
