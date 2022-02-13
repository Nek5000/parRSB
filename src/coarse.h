#ifndef _COARSE_H_
#define _COARSE_H_

#include <gslib.h>

#ifndef scalar
#define scalar double
#endif

#define CSC 0
#define CSR 1
#define IS_CSC(A) ((A)->type == CSC)
#define IS_CSR(A) ((A)->type == CSR)
#define IS_DIAG(A) ((A)->diag_val != NULL)

struct mat;
struct par_mat *par_csr_setup_ext(struct array *entries, int sep, buffer *bfr);
struct par_mat *par_csr_setup_con(const uint nelt, const ulong *eid,
                                  const slong *vtx, int nv, int sep,
                                  struct comm *c, struct crystal *cr,
                                  buffer *bfr);
void par_mat_print(struct par_mat *A);
int par_mat_free(struct par_mat *A);

struct coarse;
struct coarse *coarse_setup(uint nelt, int nv, slong const *vtx, struct comm *c,
                            buffer *bfr);
int coarse_solve(scalar *x, scalar *b, struct coarse *crs, struct comm *c,
                 buffer *bfr);
int coarse_free(struct coarse *crs);

#endif
