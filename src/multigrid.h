#ifndef _GENMAP_PRECON_H_
#define _GENMAP_PRECON_H_

#include "coarse.h"

#define TOL 1e-12

struct mg_data *mg_setup(const uint nelt, const ulong *eid, const slong *vtx,
                         int nv, int factor, struct comm *c, buffer *bfr);
void mg_vcycle(scalar *u, scalar *rhs, struct mg_data *d, struct comm *c,
               buffer *bfr);
void mg_free(struct mg_data *d);

#endif
