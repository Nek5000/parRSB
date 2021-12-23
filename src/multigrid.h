#ifndef _GENMAP_PRECON_H_
#define _GENMAP_PRECON_H_

#include <genmap-impl.h>

struct mg_lvl {
  int nsmooth;
  GenmapScalar sigma;
  struct gs_data *J; // interpolation from level i to i+1
  struct gs_data *Q; // global to local conversion of a vector
  struct csr_laplacian *M;
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

void mg_setup(struct mg_data *d, int factor, struct comm *c,
              struct csr_laplacian *M);
void mg_vcycle(GenmapScalar *u, GenmapScalar *rhs, struct mg_data *d);
void mg_free(struct mg_data *d);

#endif
