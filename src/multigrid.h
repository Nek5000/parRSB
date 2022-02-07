#ifndef _GENMAP_PRECON_H_
#define _GENMAP_PRECON_H_

#include <gslib.h>

#ifndef scalar
#define scalar double
#endif

#define TOL 1e-12

struct par_mat;
struct gs_data *setup_Q(const struct par_mat *M, const struct comm *c,
                        buffer *bfr);
void mat_vec_csr(scalar *y, scalar *x, struct par_mat *M, struct gs_data *gsh,
                 scalar *buf, buffer *bfr);

struct mg_data;
struct mg_data *mg_setup(const struct par_mat *M, const int factor,
                         const struct comm *c, struct crystal *cr, buffer *bfr);
void mg_vcycle(scalar *u, scalar *rhs, struct mg_data *d, struct comm *c,
               buffer *bfr);
void mg_free(struct mg_data *d);

#endif
