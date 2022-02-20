#ifndef _MAT_H_
#define _MAT_H_

#include <gslib.h>

#ifndef scalar
#define scalar double
#endif

struct nbr {
  ulong r, c;
  uint proc;
};

struct mat_ij {
  ulong r, c;
  uint idx, p;
  scalar v;
};

struct mat {
  ulong start, *col;
  uint n, *Lp, *Li;
  scalar *L, *D;
};

struct par_mat {
  // CSC or CSR or whatever
  int type;

  // Unique global column ids and row ids of the matrix
  uint cn, rn;
  ulong *cols, *rows;

  // Adjacency matrix
  uint *adj_off;
  uint *adj_idx;
  scalar *adj_val;

  // Diagonal
  scalar *diag_val;
};

int IS_CSC(const struct par_mat *A);
int IS_CSR(const struct par_mat *A);
int IS_DIAG(const struct par_mat *A);

void find_nbrs(struct array *arr, const ulong *eid, const slong *vtx,
               const uint nelt, const int nv, const struct comm *c,
               struct crystal *cr, buffer *buf);
int compress_nbrs(struct array *eij, struct array *nbr, buffer *bfr);

int csr_setup(struct mat *mat, struct array *entries, int sep, buffer *buf);
int mat_free(struct mat *mat);

int par_csc_setup(struct par_mat *mat, struct array *entries, int sd,
                  buffer *buf);
int par_csr_setup(struct par_mat *mat, struct array *entries, int sd,
                  buffer *buf);
struct par_mat *par_csr_setup_ext(struct array *entries, int sep, buffer *bfr);
struct par_mat *par_csr_setup_con(const uint nelt, const ulong *eid,
                                  const slong *vtx, int nv, int sep,
                                  struct comm *c, struct crystal *cr,
                                  buffer *bfr);
void par_mat_print(struct par_mat *A);
int par_mat_free(struct par_mat *A);

#endif // _MAT_H_
