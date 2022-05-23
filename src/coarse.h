#ifndef _COARSE_H_
#define _COARSE_H_

#include "gslib.h"
#include "mat.h"

struct coarse;
struct coarse *coarse_setup(const unsigned int nelt, const int nv,
                            const long long *vtx, MPI_Comm c);
int coarse_solve(scalar *x, scalar *b, struct coarse *crs, buffer *bfr);
int coarse_free(struct coarse *crs);

#endif
