#ifndef _ILU_H_
#define _ILU_H_

#include "mat.h"

struct ilu;
struct ilu *ilu_setup(const uint n, const int nv, const slong *vtx,
                      const struct comm *c, const int type, const int verbose,
                      buffer *bfr);
void ilu_free(struct ilu *ilu);

#endif
