#ifndef _COARSE_H_
#define _COARSE_H_

#include <gslib.h>

#ifndef scalar
#define scalar double
#endif

struct mat;
struct par_mat *par_csr_setup_ext(struct array *entries, int sep, buffer *bfr);
struct par_mat *par_csr_setup_con(const uint nelt, const ulong *eid,
                                  const slong *vtx, int nv, int sep,
                                  struct comm *c, struct crystal *cr,
                                  buffer *bfr);

struct coarse;
struct coarse *coarse_setup(uint nelt, int nv, slong const *vtx, struct comm *c,
                            buffer *bfr);
int coarse_solve(scalar *x, struct coarse *crs, scalar *rhs);
int coarse_free(struct coarse *crs);

#endif
