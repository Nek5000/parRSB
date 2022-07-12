#ifndef _COARSE_H_
#define _COARSE_H_

#include "gslib.h"
#include "mat.h"

struct coarse;

// API for the Laplacian (which involves solving for the dual graph)
struct coarse *coarse_setup(unsigned nelt, unsigned nv, const long long *vtx,
                            const scalar *coord, unsigned null_space,
                            unsigned type, struct comm *c);
int coarse_solve(scalar *x, scalar *b, scalar tol, struct coarse *crs,
                 buffer *bfr);
int coarse_free(struct coarse *crs);

// Alternative API for a general matrix
#define crs_schur_setup PREFIXED_NAME(crs_schur_setup)
#define crs_schur_solve PREFIXED_NAME(crs_schur_solve)
#define crs_schur_stats PREFIXED_NAME(crs_schur_stats)
#define crs_schur_free PREFIXED_NAME(crs_schur_free)

struct coarse *crs_schur_setup(uint n, const ulong *id, uint nz, const uint *Ai,
                               const uint *Aj, const scalar *A,
                               unsigned null_space, unsigned type,
                               const struct comm *comm);
void crs_schur_solve(scalar *x, struct coarse *crs, scalar *b, scalar tol);
void crs_schur_stats(struct coarse *crs);
void crs_schur_free(struct coarse *crs);

#endif
