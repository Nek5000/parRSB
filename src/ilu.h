#ifndef _ILU_H_
#define _ILU_H_

#include "mat.h"

struct ilu;
struct ilu *ilu_setup(uint n, int nv, const slong *vtx, const struct comm *c,
                      int type, buffer *bfr);
void ilu_free(struct ilu *ilu);

#endif
