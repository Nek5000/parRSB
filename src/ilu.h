#ifndef _ILU_H_
#define _ILU_H_

#include "mat.h"

struct ilu;
struct ilu *ilu_setup(const uint n, const int nv, const long long *vtx,
                      const int type, const double tol, const MPI_Comm comm,
                      const int verbose);
void ilu_free(struct ilu *ilu);

#endif
