#ifndef _COARSE_H_
#define _COARSE_H_

#include "gslib.h"
#include "mat.h"

struct coarse;
struct coarse *coarse_setup(const unsigned nelt, const unsigned nv,
                            const long long *vtx, const double *coord, int type,
                            MPI_Comm c);
int coarse_solve(scalar *x, scalar *b, scalar tol, struct coarse *crs,
                 buffer *bfr);
int coarse_free(struct coarse *crs);

// Alternative API
struct coarse *crs_setup(uint n, const ulong *id, uint nz, const uint *Ai,
                         const uint *Aj, const scalar *A, unsigned ndim,
                         scalar *coord, unsigned null_space, unsigned type,
                         const struct comm *comm);
void crs_solve(scalar *x, struct coarse *crs, scalar *b, scalar tol);
void crs_stats(struct coarse *crs);
void crs_free(struct coarse *crs);

#endif
