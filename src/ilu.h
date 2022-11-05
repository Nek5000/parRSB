#ifndef _PARRSB_ILU_H_
#define _PARRSB_ILU_H_

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

typedef struct {
  // ILU type: ILU(0), ILUC, etc.
  unsigned type;
  // Verbose level: 0, 1, etc.
  unsigned verbose;
  // Use pivoting or not: 0 or 1
  unsigned pivot;
  // 1st dropping rule: An entry a_ij is dropped abs(a_ij) < tol
  double tol;
  // 2nd dropping rule: Entries are dropped so that total nnz per row/col < p
  unsigned int nnz_per_row;
} ilu_options;

extern ilu_options ilu_default_options;

struct ilu;
struct ilu *ilu_setup(unsigned n, unsigned nv, const long long *vtx,
                      const ilu_options *options, MPI_Comm comm, buffer *bfr);
void ilu_solve(double *x, struct ilu *ilu, const double *b, buffer *bfr);
void ilu_free(struct ilu *ilu);

#ifdef __cplusplus
}
#endif

#endif
