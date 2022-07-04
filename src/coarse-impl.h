#ifndef _COARSE_IMPL_H_
#define _COARSE_IMPL_H_

#include "coarse.h"

struct coarse {
  unsigned type;       // type = schur-2-lvl, schur-3-lvl
  unsigned null_space; // Is there a null space or not
  uint un;             // User vector size
  uint cn;   // Compressed (ignoring duplicates and zero global ids) vector size
  uint an;   // Assembled size -- this is the local size of the assmebled coarse
             // matrix
  sint *u2c; // Mapping from user vector to compress vector
  struct gs_data *c2a; // Mapping from compressed vector to assmbled vector

  ulong s[3], ng[3];
  uint n[3];
  // Get rid of `idx` map -- we don't need it anymore.
  uint *idx;
  struct comm c;
  void *solver;
};

int schur_setup(struct coarse *crs, struct array *eij, struct crystal *cr,
                buffer *bfr);
int schur_solve(scalar *x, scalar *b, scalar tol, struct coarse *crs,
                buffer *bfr);
int schur_free(struct coarse *crs);

#endif
