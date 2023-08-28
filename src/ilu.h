#ifndef _PARRSB_ILU_H_
#define _PARRSB_ILU_H_

#include "mat.h"

typedef struct {
  // ILU type: ILU(0), ILUC, etc.
  unsigned type;
  // Verbose level: 0, 1, etc.
  unsigned verbose;
  // Use pivoting or not: 0 or 1
  unsigned pivot;
  // Is there a null space?
  unsigned null_space;
  // 1st dropping rule: An entry a_ij is dropped abs(a_ij) < tol
  scalar tol;
  // 2nd dropping rule: Entries are dropped so that total nnz per row/col < p
  unsigned int nnz_per_row;
} ilu_options;

struct ilu;
struct ilu *ilu_setup(unsigned n, unsigned nv, const long long *vtx,
                      const ilu_options *options, MPI_Comm comm, buffer *bfr);
void ilu_solve(double *x, const struct ilu *ilu, const double *b, buffer *bfr);
void ilu_free(struct ilu *ilu);

#endif
