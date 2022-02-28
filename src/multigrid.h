#ifndef _MULTIGRID_H_
#define _MULTIGRID_H_

#include "mat.h"

struct mg_data;
struct mg_data *mg_setup(const struct par_mat *M, const int factor,
                         const struct comm *c, struct crystal *cr, buffer *bfr);
void mg_vcycle(scalar *u, scalar *rhs, struct mg_data *d, struct comm *c,
               buffer *bfr);
void mg_free(struct mg_data *d);

#endif
