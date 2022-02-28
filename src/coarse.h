#ifndef _COARSE_H_
#define _COARSE_H_

#include "gslib.h"
#include "mat.h"

struct coarse {
  ulong ls, lg, is, ig;
  uint ln, in, *idx;
  void *solver;
};

struct coarse *coarse_setup(uint nelt, int nv, slong const *vtx, struct comm *c,
                            buffer *bfr);
int coarse_solve(scalar *x, scalar *b, struct coarse *crs, struct comm *c,
                 buffer *bfr);
int coarse_free(struct coarse *crs);

#endif
