#ifndef _PARRSB_COARSE_H_
#define _PARRSB_COARSE_H_

#include "gslib.h"

#if !defined(MPI)
#error "gslib needs to be compiled with MPI"
#endif

#if !defined(GLOBAL_LONG_LONG)
#error "gslib needs to be compiled with GLOBAL_LONG_LONG"
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct coarse;

// API for the Laplacian (which involves solving for the dual graph)
struct coarse *coarse_setup(unsigned nelt, unsigned nv, const long long *vtx,
                            const double *coord, unsigned null_space,
                            unsigned type, struct comm *c);
void coarse_solve(double *x, struct coarse *crs, double *b, double tol);
void coarse_free(struct coarse *crs);

// Alternative API for a general matrix
#define crs_parrsb_setup PREFIXED_NAME(crs_parrsb_setup)
struct coarse *crs_parrsb_setup(uint n, const ulong *id, uint nz,
                                const uint *Ai, const uint *Aj, const double *A,
                                unsigned null_space, unsigned type,
                                const struct comm *comm);

#define crs_parrsb_solve PREFIXED_NAME(crs_parrsb_solve)
void crs_parrsb_solve(double *x, struct coarse *crs, double *b, double tol);

#define crs_parrsb_free PREFIXED_NAME(crs_parrsb_free)
void crs_parrsb_free(struct coarse *crs);

#ifdef __cplusplus
}
#endif

#endif
