#ifndef _COARSE_H_
#define _COARSE_H_

#include <gslib.h>

#define scalar double

struct coarse;
struct coarse *coarse_setup(uint nelt, int nv, slong const *vtx, struct comm *c,
                            buffer *bfr);
int coarse_solve(scalar *x, struct coarse *crs, scalar *rhs);
int coarse_free(struct coarse *crs);

#endif
