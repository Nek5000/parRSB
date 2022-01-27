#ifndef _COARSE_H_
#define _COARSE_H_

#include <gslib.h>

#define scalar double

struct coarse;
struct coarse *coarse_setup(uint nelt, int nv, slong const *vtx, struct comm *c,
                            buffer *bfr);
int coarse_solve(scalar *x, struct coarse *crs, scalar *rhs);
int coarse_free(struct coarse *crs);

struct mat;
struct mat *csr_setup_ext(struct array *entries, int sep, buffer *bfr);
struct mat *csr_from_conn(const uint nelt, const ulong *eid, const slong *vtx,
                          int nv, int sep, struct comm *c, buffer *bfr);
int mat_zero_sum(struct mat *mat);
int mat_print(struct mat *mat);
int mat_free(struct mat *mat);

#endif
