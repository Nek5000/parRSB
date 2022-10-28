#ifndef _PARRSB_MAT_H_
#define _PARRSB_MAT_H_

#include "gslib.h"

#ifdef scalar
#undef scalar
#endif
#define scalar double

#ifdef SCALAR_MAX
#undef SCALAR_MAX
#endif
#define SCALAR_MAX DBL_MAX

//==============================================================================
// Memory management
//
int sfree(void *p, const char *file, unsigned line);
#define tfree(x) sfree(x, __FILE__, __LINE__)

//==============================================================================
// Helper functions
//
struct nbr {
  ulong r, c;
  uint proc;
};

// `arr` is an array of type `struct nbr`
void find_nbrs(struct array *arr, const ulong *eid, const slong *vtx,
               const uint nelt, const int nv, struct crystal *cr, buffer *buf);

struct mij {
  ulong r, c;
  uint idx, p;
  scalar v;
};

// `eij` is an array of type `struct mij`, input `nbr` is an array of type
// `struct nbr`
int compress_nbrs(struct array *eij, struct array *nbr, buffer *bfr);

//==============================================================================
// mat
// * `entries` is an array of type `struct mij`
// * sd = 0 (don't store diagonal separately) or 1 (store diagonal separately)

struct mat {
  ulong start;
  uint n, *Lp, *Li;
  scalar *L, *D;
};

int mat_setup(struct mat *A, struct array *entries, int sd, buffer *buf);
void mat_print(const struct mat *A);
void mat_dump(const char *name, const struct mat *A, struct crystal *cr,
              buffer *bfr);
int mat_free(struct mat *A);

//==============================================================================
// par_mat
// * entries is an array of type `struct mij`
// * type = 0 (CSC) or 1 (CSR)
// * sd = 0 (don't store diagonal separately) or 1 (store diagonal separately)

struct par_mat {
  int type;
  uint cn, rn, *adj_off, *adj_idx, *diag_idx;
  scalar *adj_val, *diag_val;
  ulong *cols, *rows;
};

int IS_CSC(const struct par_mat *A);
int IS_CSR(const struct par_mat *A);
int IS_DIAG(const struct par_mat *A);

int par_mat_setup(struct par_mat *M, struct array *entries, int type, int sd,
                  buffer *bfr);
int par_csc_setup(struct par_mat *mat, struct array *entries, int sd,
                  buffer *buf);
int par_csr_setup(struct par_mat *mat, struct array *entries, int sd,
                  buffer *buf);
void par_csr_to_csc(struct par_mat *B, const struct par_mat *A, int sd,
                    struct crystal *cr, buffer *bfr);
void par_csc_to_csr(struct par_mat *B, const struct par_mat *A, int sd,
                    struct crystal *cr, buffer *bfr);
void par_mat_print(struct par_mat *A);
void par_mat_dump(const char *name, struct par_mat *A, struct crystal *cr,
                  buffer *bfr);
int par_mat_free(struct par_mat *A);

// Create a par_mat from connectivity
struct par_mat *par_csr_setup_con(const uint nelt, const ulong *eid,
                                  const slong *vtx, int nv, int sep,
                                  struct comm *c, struct crystal *cr,
                                  buffer *bfr);
struct par_mat *par_csr_setup_ext(struct array *entries, int sd, buffer *bfr);

// Mat vec routines
struct gs_data *setup_Q(const struct par_mat *M, const struct comm *c,
                        buffer *bfr);
int par_mat_vec(scalar *y, const scalar *x, const struct par_mat *M,
                struct gs_data *gsh, scalar *buf, buffer *bfr);

#endif // _PARRSB_MAT_H_
