#include "coarse.h"
#include "multigrid.h"
#include <math.h>

//------------------------------------------------------------------------------
// Cholesky factorization of a mat
//
/*
symbolic factorization: finds the sparsity structure of L

uses the concept of elimination tree:
  the parent of node j is node i when L(i,j) is the first
    non-zero in column j below the diagonal (i>j)
  L's structure is discovered row-by-row; the first time
    an entry in column j is set, it must be the parent

the nonzeros in L are the nonzeros in A + paths up the elimination tree

linear in the number of nonzeros of L
*/
static uint *cholesky_symbolic(struct mat *L, uint n, uint const *Ap,
                               uint const *Ai, buffer *buf) {
  L->n = n;

  uint *parent = tcalloc(uint, 2 * n), *visit = parent + n;
  uint i, j, nz = 0;
  for (i = 0; i < n; i++) {
    parent[i] = n, visit[i] = i;
    for (uint p = Ap[i]; p < Ap[i + 1]; p++) {
      if ((j = Ai[p]) >= i)
        break;
      for (; visit[j] != i; j = parent[j]) {
        ++nz, visit[j] = i;
        if (parent[j] == n) {
          parent[j] = i;
          break;
        }
      }
    }
  }

  uint *Lp = L->Lp = tcalloc(uint, n + 1);
  uint *Li = L->Li = tcalloc(uint, nz);

  Lp[0] = 0;
  uint *Lir, nzr;
  for (i = 0; i < n; i++) {
    nzr = 0, Lir = &Li[Lp[i]];
    visit[i] = i;
    for (uint p = Ap[i]; p < Ap[i + 1]; p++) {
      if ((j = Ai[p]) >= i)
        break;
      for (; visit[j] != i; j = parent[j])
        Lir[nzr++] = j, visit[j] = i;
    }
    sortv(Lir, Lir, nzr, sizeof(uint), buf);
    Lp[i + 1] = Lp[i] + nzr;
  }

  free(parent);
}

/*
numeric factorization:

L is built row-by-row, using:    ( ' indicates transpose )


[ A  r ]  = [ (I-L)   ] [ D^(-1)  ] [ (I-L)' -s ]
[ r' a ]    [  -s'  1 ] [     1/d ] [         1 ]

          = [ A   (I-L) D^(-1) (-s)  ]
            [ r'  s' D^(-1) s + 1/d  ]

so, if r' is the next row of A, up to but excluding the diagonal,
then the next row of L, s', obeys

   r = - (I-L) D^(-1) s

let y = (I-L)^(-1) (-r)
then s = D y, and d = 1/(a - s' y)
*/
static void cholesky_numeric(struct mat *chol, const uint n, const uint *Ap,
                             const uint *Ai, const scalar *A, uint *visit,
                             scalar *y) {
  const uint *Lp = chol->Lp, *Li = chol->Li;
  scalar *D = chol->D = tcalloc(scalar, n);
  scalar *L = chol->L = tcalloc(scalar, Lp[n]);

  uint i;
  for (i = 0; i < n; i++) {
    uint p, pe, j;
    scalar a;
    visit[i] = n;
    for (p = Lp[i], pe = Lp[i + 1]; p != pe; p++)
      j = Li[p], y[j] = 0, visit[j] = i;
    for (p = Ap[i], pe = Ap[i + 1]; p != pe; p++) {
      if ((j = Ai[p]) >= i) {
        if (j == i)
          a = A[p];
        break;
      }
      y[j] = -A[p];
    }
    for (p = Lp[i], pe = Lp[i + 1]; p != pe; p++) {
      uint j = Li[p], q = Lp[j], qe = Lp[j + 1];
      scalar yj = y[j];
      for (; q != qe; q++) {
        uint k = Li[q];
        if (visit[k] == i)
          yj += L[q] * y[k];
      }
      y[j] = yj;
      scalar lij = L[p] = D[j] * yj;
      a -= lij * yj;
    }
    D[i] = 1 / a;
  }
}

static void cholesky_factor(struct mat *L, struct mat *A, buffer *buf) {
  L->start = A->start;
  const uint uints_as_dbls =
      (A->n * sizeof(uint) + sizeof(double) - 1) / sizeof(double);
  buffer_reserve(buf, (uints_as_dbls + A->n) * sizeof(double));
  cholesky_symbolic(L, A->n, A->Lp, A->Li, buf);
  cholesky_numeric(L, L->n, A->Lp, A->Li, A->L, buf->ptr,
                   uints_as_dbls + (double *)buf->ptr);
}

static void cholesky_solve(scalar *x, const struct mat *A, scalar *b) {
  const uint *Lp = A->Lp, *Li = A->Li, n = A->n;
  const scalar *L = A->L, *D = A->D;

  uint i, p, pe;
  for (i = 0; i < n; i++) {
    scalar xi = b[i];
    for (p = Lp[i], pe = Lp[i + 1]; p != pe; p++)
      xi += x[Li[p]] * L[p];
    x[i] = xi;
  }

  for (i = 0; i < n; i++)
    x[i] *= D[i];

  for (i = n; i > 0;) {
    scalar xi = x[--i];
    for (p = Lp[i], pe = Lp[i + 1]; p != pe; p++) {
      x[Li[p]] += xi * L[p];
    }
    x[i] = xi;
  }
}

static void cholesky_lower_solve(scalar *x, const struct mat *A, scalar *b) {
  const uint *Lp = A->Lp, *Li = A->Li, n = A->n;
  const scalar *L = A->L, *D = A->D;

  uint i, p, pe;
  for (i = 0; i < n; i++) {
    scalar xi = b[i];
    for (p = Lp[i], pe = Lp[i + 1]; p != pe; p++)
      xi += x[Li[p]] * L[p];
    x[i] = xi;
  }

  for (i = 0; i < n; i++)
    x[i] *= sqrt(D[i]);
}

static void cholesky_upper_solve(scalar *x, const struct mat *A, scalar *b) {
  const uint *Lp = A->Lp, *Li = A->Li, n = A->n;
  const scalar *L = A->L, *D = A->D;

  uint i;
  for (i = 0; i < n; i++)
    x[i] = b[i] * sqrt(D[i]);

  uint p, pe;
  for (i = n; i > 0;) {
    scalar xi = x[--i];
    for (p = Lp[i], pe = Lp[i + 1]; p != pe; p++) {
      x[Li[p]] += xi * L[p];
    }
    x[i] = xi;
  }
}

//-----------------------------------------------------------------------------
// Schur setup, solve and free
//
// A_ll: local dof of a processor (block diagonal across processors)
// A_sl: shared - local matrix
// A_ss: shared dof freedom (matrix is split row wise)
//
//     |A_ll (B)  A_ls (F)|
//  A= |                  |
//     |A_sl (E)  A_ss (S)|
//
struct schur {
  struct mat A_ll;
  struct par_mat A_ls, A_sl, A_ss;
  struct gs_data *Q_ls, *Q_sl, *Q_ss;
  struct mg *M;
};

static int S_owns_row(const ulong r, const ulong *rows, const uint n) {
  // We can do a binary search instead of linear search
  uint i = 0;
  while (i < n && rows[i] != r)
    i++;
  return i;
}

// Calculate G = L_{B}^{-1} x F where B = L_{B} U_{B}. F is in CSR format,
// distributed by rows similar to B. G will be in CSC format and distributed
// by columns similar to row distribution of S.
static int schur_setup_G(struct par_mat *G, scalar tol, const struct mat *L,
                         const struct par_mat *F, const ulong *srows,
                         const uint srn, const struct comm *c,
                         struct crystal *const cr, buffer *bfr) {
  assert(IS_CSR(F));
  assert(!IS_DIAG(F));

  buffer_reserve(bfr, sizeof(scalar) * L->n * F->cn);
  scalar *v = (scalar *)bfr->ptr;
  for (uint i = 0; i < L->n * F->cn; i++)
    v[i] = 0;

  scalar *b = tcalloc(scalar, 2 * L->n);
  scalar *x = b + L->n;
  for (uint i = 0; i < F->rn; i++) {
    b[F->rows[i] - L->start] = 1;
    cholesky_lower_solve(x, L, b);

    // Calculate F: i^th row of F is multiplied by each element of i^th
    // column of L_B^-1
    for (uint k = F->adj_off[i], ke = F->adj_off[i + 1]; k < ke; k++)
      for (uint j = 0; j < L->n; j++)
        // m.c = F->cols[F->adj_idx[k]], m.r = L->start + j;
        v[j * F->cn + F->adj_idx[k]] += F->adj_val[k] * x[j];

    b[F->rows[i] - L->start] = 0;
    for (uint j = 0; j < L->n; j++)
      x[j] = 0;
  }

  uint size = L->n * 20 + 1, k;
  struct mij m;

  struct array unique;
  array_init(struct mij, &unique, size);

  slong *cols = tcalloc(slong, size);
  sint *owners = tcalloc(sint, size);
  for (uint i = k = 0; i < L->n; i++) {
    for (uint j = 0; j < F->cn; j++) {
      if (fabs(v[i * F->cn + j]) >= tol) {
        m.r = L->start + i, m.c = F->cols[j], m.v = v[i * F->cn + j];
        if (k == size) {
          size += size / 2 + 1;
          cols = (slong *)trealloc(slong, cols, size);
          owners = (sint *)trealloc(sint, owners, size);
        }
        owners[k] = S_owns_row(m.c, srows, srn) < srn ? c->id : -1;
        cols[k++] = m.c;
        array_cat(struct mij, &unique, &m, 1);
      }
    }
  }

  // Set the destination processor for elements in unique. Since W will be in
  // CSR and share the same row distribution as S, we use gslib and row
  // distribution of S to find the destination processor.
  struct gs_data *gsh = gs_setup(cols, k, c, 0, gs_pairwise, 1);
  gs(owners, gs_int, gs_max, 0, gsh, bfr);
  free(gsh);

  struct mij *ptr = unique.ptr;
  for (uint i = 0; i < unique.n; i++)
    ptr[i].p = owners[i];
  free(owners), free(cols);

  sarray_transfer(struct mij, &unique, p, 0, cr);

  par_csc_setup(G, &unique, 0, bfr);
#ifdef DUMPG
  par_mat_print(G);
#endif
  array_free(&unique);

  return 0;
}

// Calculate W = E x U_{B}^{-1} where B = L_{B} U_{B}. E is in CSC format.
// W will be in CSR format and distributed by rows similar to distribution of S.
static int schur_setup_W(struct par_mat *W, scalar tol, const struct mat *L,
                         const struct par_mat *E, const ulong *srows,
                         const uint srn, const struct comm *c,
                         struct crystal *const cr, buffer *bfr) {
  assert(IS_CSC(E));
  assert(!IS_DIAG(E));

  buffer_reserve(bfr, sizeof(scalar) * L->n * E->rn);
  scalar *v = (scalar *)bfr->ptr;
  for (uint i = 0; i < L->n * E->rn; i++)
    v[i] = 0;

  // Multiply E by U_B^{-1} now. Columns of U_B^{-1} are found one by one and
  // then E is multiplied by each column.
  scalar *b = tcalloc(scalar, 2 * L->n);
  scalar *x = b + L->n;
  for (uint i = 0; i < L->n; i++) {
    b[i] = 1;
    cholesky_upper_solve(x, L, b);

    // Multiply E by x: i^th col of E is multiplied by element x[i]
    for (uint j = 0; j < E->cn; j++)
      for (uint k = E->adj_off[j], ke = E->adj_off[j + 1]; k < ke; k++)
        // m.c = L->start + i, m.r = E->rows[E->adj_idx[k]];
        v[E->adj_idx[k] * L->n + i] += E->adj_val[k] * x[E->cols[j] - L->start];

    b[i] = 0;
    for (uint j = 0; j < L->n; j++)
      x[j] = 0;
  }

  uint size = E->rn * 20 + 1, k;
  struct mij m;

  struct array unique;
  array_init(struct mij, &unique, size);

  slong *rows = (slong *)tcalloc(slong, size);
  sint *owners = (sint *)tcalloc(sint, size);
  for (uint i = k = 0; i < E->rn; i++) {
    for (uint j = 0; j < L->n; j++) {
      if (fabs(v[i * L->n + j]) > tol) {
        m.r = E->rows[i], m.c = L->start + j, m.v = v[i * L->n + j];
        if (k == size) {
          size += size / 2 + 1;
          rows = (slong *)trealloc(slong, rows, size);
          owners = (sint *)trealloc(sint, owners, size);
        }
        owners[k] = S_owns_row(m.r, srows, srn) < srn ? c->id : -1;
        rows[k++] = m.r;
        array_cat(struct mij, &unique, &m, 1);
      }
    }
  }

  // Set the destination processor for elements in unique. Since W will be in
  // CSR and share the same row distribution as S, we use gslib and row
  // distribution of S to find the destination processor.
  struct gs_data *gsh = gs_setup(rows, k, c, 0, gs_pairwise, 0);
  gs(owners, gs_int, gs_max, 0, gsh, bfr);
  free(gsh);

  struct mij *ptr = unique.ptr;
  for (uint i = 0; i < unique.n; i++)
    ptr[i].p = owners[i];
  free(owners), free(rows);

  sarray_transfer(struct mij, &unique, p, 0, cr);

  par_csr_setup(W, &unique, 0, bfr);
#ifdef DUMPW
  par_mat_print(W);
#endif
  array_free(&unique);

  return 0;
}

// C = A - B; A and B should be in CSR format with the same row
// distribution across processors
static int sparse_sub(struct par_mat *C, const struct par_mat *A,
                      const struct par_mat *B, const struct comm *c,
                      struct crystal *const cr, buffer *bfr) {
  assert(IS_CSR(A));
  assert(IS_CSR(B));

  struct array cij;
  array_init(struct mij, &cij, 100);

  struct mij m;
  uint r, j, je;
  for (r = 0; r < A->rn; r++) {
    m.r = A->rows[r];
    for (j = A->adj_off[r], je = A->adj_off[r + 1]; j != je; j++) {
      m.c = A->cols[A->adj_idx[j]], m.v = A->adj_val[j];
      array_cat(struct mij, &cij, &m, 1);
    }
  }
  if (IS_DIAG(A)) {
    for (r = 0; r < A->rn; r++) {
      m.r = A->rows[r], m.c = A->rows[r], m.v = A->diag_val[r];
      array_cat(struct mij, &cij, &m, 1);
    }
  }

  for (r = 0; r < B->rn; r++) {
    m.r = B->rows[r];
    for (j = B->adj_off[r], je = B->adj_off[r + 1]; j != je; j++) {
      m.c = B->cols[B->adj_idx[j]], m.v = -B->adj_val[j];
      array_cat(struct mij, &cij, &m, 1);
    }
  }
  if (IS_DIAG(B)) {
    for (r = 0; r < B->rn; r++) {
      m.r = B->rows[r], m.c = B->rows[r], m.v = -B->diag_val[r];
      array_cat(struct mij, &cij, &m, 1);
    }
  }

  sarray_sort_2(struct mij, cij.ptr, cij.n, r, 1, c, 1, bfr);
  struct mij *ptr = (struct mij *)cij.ptr;

  struct array unique;
  array_init(struct mij, &unique, 100);

  if (cij.n > 0) {
    uint i = 0;
    while (i < cij.n) {
      scalar s = 0;
      for (j = i; j < cij.n && ptr[j].r == ptr[i].r && ptr[j].c == ptr[i].c;
           j++)
        s += ptr[j].v;
      m = ptr[i], m.v = s;
      array_cat(struct mij, &unique, &m, 1);
      i = j;
    }
  }
  array_free(&cij);

  par_csr_setup(C, &unique, 1, bfr);
  array_free(&unique);

  return 0;
}

static int sparse_gemm(struct par_mat *WG, const struct par_mat *W,
                       const struct par_mat *G, struct crystal *cr,
                       const struct comm *c, buffer *bfr) {
  // W is in CSR, G is in CSC; we multiply rows of W by shifting
  // the columns of G from processor to processor. This is not scalable
  // at all -- need to do a 2D partition of the matrices W and G.
  assert(IS_CSR(W));
  assert(IS_CSC(G));

  // Put G into an array to transfer from processor to processor
  struct array gij, sij;
  array_init(struct mij, &gij, 100);
  array_init(struct mij, &sij, 100);

  struct mij m;
  uint i, j, je;
  for (i = 0; i < G->cn; i++) {
    m.c = G->cols[i], m.p = c->id;
    for (j = G->adj_off[i], je = G->adj_off[i + 1]; j != je; j++) {
      m.r = G->rows[G->adj_idx[j]];
      m.v = G->adj_val[j];
      array_cat(struct mij, &gij, &m, 1);
    }
  }

  sarray_sort_2(struct mij, gij.ptr, gij.n, r, 1, c, 1, bfr);
  struct mij *ptr = (struct mij *)gij.ptr;
  for (i = 0; i < gij.n; i++)
    ptr[i].idx = i;

  uint idx;
  for (uint p = 0; p < c->np; p++) {
    for (i = 0; i < W->rn; i++) {
      m.r = W->rows[i], idx = 0;
      for (j = W->adj_off[i], je = W->adj_off[i + 1]; j < je; j++) {
        ulong k = W->cols[W->adj_idx[j]];
        while (idx < gij.n && ptr[idx].r < k)
          idx++;
        while (idx < gij.n && ptr[idx].r == k) {
          m.c = ptr[idx].c;
          m.v = W->adj_val[j] * ptr[idx].v, idx++;
          array_cat(struct mij, &sij, &m, 1);
        }
      }
    }

    sint next = (c->id + 1) % c->np;
    for (i = 0; i < gij.n; i++)
      ptr[i].p = next;
    sarray_transfer(struct mij, &gij, p, 0, cr);

    sarray_sort(struct mij, gij.ptr, gij.n, idx, 0, bfr);
    ptr = gij.ptr;
  }

  par_csr_setup(WG, &sij, 0, bfr);

  array_free(&gij);
  array_free(&sij);

  return 0;
}

static struct mg *
schur_precond_setup(const struct mat *L, const struct par_mat *F,
                    const struct par_mat *S, const struct par_mat *E,
                    const struct comm *c, struct crystal *cr, buffer *bfr) {
  // TODO: Sparsify W and G when they are built
  struct par_mat W, G, WG;
  schur_setup_G(&G, 0.1, L, F, S->rows, S->rn, c, cr, bfr);
  schur_setup_W(&W, 0.1, L, E, S->rows, S->rn, c, cr, bfr);
  sparse_gemm(&WG, &W, &G, cr, c, bfr);
#ifdef DUMPWG
  par_mat_print(&WG);
#endif

  // P is CSR
  struct par_mat *P = tcalloc(struct par_mat, 1);
  sparse_sub(P, S, &WG, c, cr, bfr);
#ifdef DUMPP
  par_mat_print(P);
#endif

  par_mat_free(&W);
  par_mat_free(&G);
  par_mat_free(&WG);

  return mg_setup(P, 2, c, cr, bfr);
}

static struct gs_data *setup_Ezl_Q(struct par_mat *E, ulong s, uint n,
                                   struct comm *c, buffer *bfr) {
  assert(IS_CSC(E));
  assert(!IS_DIAG(E));

  buffer_reserve(bfr, sizeof(slong) * (n + E->rn));
  slong *ids = (slong *)bfr->ptr;
  uint i, j;
  for (i = 0; i < n; i++)
    ids[i] = s + i;
  for (j = 0; j < E->rn; j++, i++)
    ids[i] = -E->rows[j];

  return gs_setup(ids, i, c, 0, gs_pairwise, 0);
}

static int Ezl(scalar *y, const struct par_mat *E, struct gs_data *gsh,
               scalar *zl, const ulong s, const uint n, buffer *bfr) {
  assert(IS_CSC(E));
  assert(!IS_DIAG(E));

  uint nn = n + E->rn;
  scalar *wrk = (scalar *)tcalloc(scalar, nn);
  scalar *ye = wrk + n;
  for (uint i = 0; i < E->cn; i++) {
    scalar zlk = zl[E->cols[i] - s];
    for (uint j = E->adj_off[i], je = E->adj_off[i + 1]; j < je; j++)
      ye[E->adj_idx[j]] += zlk * E->adj_val[j];
  }

  gs(wrk, gs_double, gs_add, 1, gsh, bfr);

  for (uint i = 0; i < n; i++)
    y[i] = wrk[i];

  free(wrk);

  return 0;
}

static struct gs_data *setup_Fxi_Q(struct par_mat *F, ulong s, uint n,
                                   struct comm *c, buffer *bfr) {
  assert(IS_CSR(F));
  assert(!IS_DIAG(F));

  uint nnz = F->rn > 0 ? F->adj_off[F->rn] : 0;
  buffer_reserve(bfr, sizeof(slong) * (n + nnz));
  slong *ids = (slong *)bfr->ptr;
  uint i, j;
  for (i = 0; i < nnz; i++)
    ids[i] = F->cols[F->adj_idx[i]];
  for (j = 0; j < n; j++, i++)
    ids[i] = -(s + j);

  return gs_setup(ids, i, c, 0, gs_pairwise, 0);
}

static int Fxi(scalar *y, const struct par_mat *F, struct gs_data *gsh,
               scalar *xi, const ulong s, const uint n, buffer *bfr) {
  assert(IS_CSR(F));
  assert(!IS_DIAG(F));

  uint nnz = F->rn > 0 ? F->adj_off[F->rn] : 0;
  scalar *wrk = (scalar *)tcalloc(scalar, nnz + n);
  uint i, j;
  for (i = 0; i < nnz; i++)
    wrk[i] = 0;
  for (j = 0; j < n; j++, i++)
    wrk[i] = xi[j];

  gs(wrk, gs_double, gs_add, 1, gsh, bfr);

  for (i = 0; i < F->rn; i++) {
    scalar si = 0;
    for (uint j = F->adj_off[i], je = F->adj_off[i + 1]; j < je; j++)
      si += F->adj_val[j] * wrk[j];
    y[F->rows[i] - s] = si;
  }
}

static int distribute_by_columns(struct array *aij, ulong s, uint n,
                                 struct comm *c, struct crystal *cr,
                                 buffer *bfr) {
  buffer_reserve(bfr, sizeof(slong) * aij->n);
  slong *cols = (slong *)bfr->ptr;
  sint *owner = (sint *)tcalloc(sint, aij->n);

  uint i;
  struct mij *ptr = aij->ptr;
  for (i = 0; i < aij->n; i++) {
    cols[i] = ptr[i].c;
    if (ptr[i].c >= s && ptr[i].c < s + n)
      owner[i] = c->id;
    else
      owner[i] = -1;
  }

  struct gs_data *gsh = gs_setup(cols, aij->n, c, 0, gs_pairwise, 0);
  gs(owner, gs_int, gs_max, 0, gsh, bfr);

  for (i = 0; i < aij->n; i++)
    ptr[i].p = owner[i];
  free(owner);

  sarray_transfer(struct mij, aij, p, 1, cr);

  return 0;
}

static inline scalar dot(scalar *r, scalar *s, uint n) {
  scalar t = 0;
  for (uint i = 0; i < n; i++)
    t += r[i] * s[i];
  return t;
}

static inline void ortho(scalar *q, uint n, ulong ng, struct comm *c) {
  scalar s = 0, buf;
  for (uint i = 0; i < n; i++)
    s += q[i];

  comm_allreduce(c, gs_double, gs_add, &s, 1, &buf);
  s /= ng;

  for (uint i = 0; i < n; i++)
    q[i] -= s;
}

static int schur_action(scalar *y, const struct schur *schur, scalar *x,
                        ulong ls, scalar *wrk, buffer *bfr) {
  const struct par_mat *S = &schur->A_ss;
  assert(IS_CSR(S));
  assert(S->rn == 0 || IS_DIAG(S));

  uint ln = schur->A_ll.n, in = S->rn;
  uint mn = ln > in ? ln : in;
  scalar *xl = (scalar *)tcalloc(scalar, 2 * mn), *exl = xl + mn;

  // Calculate (E (B^-1) F) x
  // Fx: x has size in, Fx has size ln. So wrk has to be at least ln
  Fxi(exl, &schur->A_ls, schur->Q_ls, x, ls, in, bfr);
  // Multiply Fx by B^-1 or (LU)^-1
  cholesky_solve(xl, &schur->A_ll, exl);
  // Multuply (B^-1)Fx by E
  Ezl(exl, &schur->A_sl, schur->Q_sl, xl, ls, in, bfr);

  // Separately calculate Sx
  mat_vec_csr(y, x, S, schur->Q_ss, wrk, bfr);

  for (uint i = 0; i < in; i++)
    y[i] -= exl[i];

  free(xl);

  return 0;
}

static int project(scalar *x, scalar *b, const struct schur *schur, ulong ls,
                   struct comm *c, int miter, int null_space, int verbose,
                   buffer *bfr) {
  const struct par_mat *S = &schur->A_ss;
  assert(IS_CSR(S));
  assert(S->rn == 0 || IS_DIAG(S));

  struct mg *d = schur->M;

  slong out[2][1], buf[2][1], in = S->rn;
  comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
  ulong ng = out[1][0];

  if (ng == 0)
    return 0;

  uint n = S->rn, nnz = n > 0 ? S->adj_off[n] + n : 0;
  scalar *z = (scalar *)tcalloc(scalar, 6 * n + nnz);
  scalar *w = z + n, *r = w + n, *p = r + n, *z0 = p + n, *dz = z0 + n;
  scalar *wrk = dz + n;
  scalar *P = (scalar *)tcalloc(scalar, 2 * (miter + 1) * n);
  scalar *W = P + n * (miter + 1);

  uint i;
  for (i = 0; i < n; i++)
    x[i] = 0, r[i] = b[i];

  scalar rr = dot(r, r, n);
  comm_allreduce(c, gs_double, gs_add, &rr, 1, buf);
  scalar tol = 1e-3, rtol = rr * tol * tol;

  for (i = 0; i < n; i++)
    z[i] = r[i];
  if (null_space)
    ortho(z, n, ng, c);
  scalar rz1 = dot(z, z, n);
  comm_allreduce(c, gs_double, gs_add, &rz1, 1, buf);

  for (i = 0; i < n; i++)
    p[i] = z[i];

  scalar alpha, beta, rzt, rz2;

  uint j, k;
  for (i = 0; i < miter; i++) {
    // Action of S - E (LU)^-1 F
    schur_action(w, schur, p, ls, wrk, bfr);

    scalar pw = dot(p, w, n);
    comm_allreduce(c, gs_double, gs_add, &pw, 1, buf);
    alpha = rz1 / pw;

    pw = 1 / sqrt(pw);
    for (j = 0; j < n; j++)
      W[i * n + j] = pw * w[j], P[i * n + j] = pw * p[j];

    for (j = 0; j < n; j++)
      x[j] += alpha * p[j], r[j] -= alpha * w[j];

    rr = dot(r, r, n);
    comm_allreduce(c, gs_double, gs_add, &rr, 1, buf);

    if (rr < rtol || sqrt(rr) < tol)
      break;

    for (j = 0; j < n; j++)
      z0[j] = z[j];

#if 1
    mg_vcycle(z, r, d, c, bfr);
#else
    for (j = 0; j < n; j++)
      z[j] = r[j];
#endif

    rzt = rz1;
    if (null_space)
      ortho(z, n, ng, c);
    rz1 = dot(r, z, n);
    comm_allreduce(c, gs_double, gs_add, &rz1, 1, buf);

    for (j = 0; j < n; j++)
      dz[j] = z[j] - z0[j];
    rz2 = dot(r, dz, n);
    comm_allreduce(c, gs_double, gs_add, &rz2, 1, buf);

    if (c->id == 0 && verbose > 0)
      printf("rr = %e rtol = %e rz0 = %e rz1 = %e rz2 = %e\n", rr, rtol, rzt,
             rz1, rz2);

    beta = rz2 / rzt;
    for (j = 0; j < n; j++)
      p[j] = z[j] + beta * p[j];

    for (k = 0; k < n; k++)
      P[miter * n + k] = 0;

    for (j = 0; j <= i; j++) {
      pw = 0;
      for (k = 0; k < n; k++)
        pw += W[j * n + k] * p[k];
      comm_allreduce(c, gs_double, gs_add, &pw, 1, buf);
      for (k = 0; k < n; k++)
        P[miter * n + k] += pw * P[j * n + k];
    }

    for (k = 0; k < n; k++)
      p[k] -= P[miter * n + k];
  }

  free(z);
  free(P);

  return i == miter ? i : i + 1;
}

int schur_setup(struct coarse *crs, struct array *eij, const ulong ng,
                struct comm *c, struct crystal *cr, buffer *bfr) {
  struct schur *schur = crs->solver = (struct schur *)tcalloc(struct schur, 1);

  // Setup A_ll
  struct array ll, ls, sl, ss;
  array_init(struct mij, &ll, eij->n / 4 + 1);
  array_init(struct mij, &ls, eij->n / 4 + 1);
  array_init(struct mij, &sl, eij->n / 4 + 1);
  array_init(struct mij, &ss, eij->n / 4 + 1);

  struct mij *ptr = (struct mij *)eij->ptr;
  for (uint i = 0; i < eij->n; i++) {
    if (ptr[i].r <= ng) {
      if (ptr[i].c <= ng)
        array_cat(struct mij, &ll, &ptr[i], 1);
      else
        array_cat(struct mij, &ls, &ptr[i], 1);
    } else if (ptr[i].c <= ng)
      array_cat(struct mij, &sl, &ptr[i], 1);
    else
      array_cat(struct mij, &ss, &ptr[i], 1);
  }
  assert(ls.n == sl.n); // For symmetric matrices

  // Setup local block diagonal (B)
  struct mat B;
  csr_setup(&B, &ll, 0, bfr);
  cholesky_factor(&schur->A_ll, &B, bfr);
#ifdef DUMPA
  mat_print(&B);
#endif
  mat_free(&B);
  array_free(&ll);
  schur->A_ll.start = crs->ls;

  // Setup E
  assert(schur->A_ll.n == crs->ln);
  distribute_by_columns(&sl, crs->ls, crs->ln, c, cr, bfr);
  par_csc_setup(&schur->A_sl, &sl, 0, bfr);
#ifdef DUMPA
  par_mat_print(&schur->A_sl);
#endif
  array_free(&sl);
  schur->Q_sl = setup_Ezl_Q(&schur->A_sl, crs->lg + crs->is, crs->in, c, bfr);

  // Setup S
  par_csr_setup(&schur->A_ss, &ss, 1, bfr);
#ifdef DUMPA
  par_mat_print(&schur->A_ss);
#endif
  array_free(&ss);
  schur->Q_ss = setup_Q(&schur->A_ss, c, bfr);

  // Setup F
  par_csr_setup(&schur->A_ls, &ls, 0, bfr);
#ifdef DUMPA
  par_mat_print(&schur->A_ls);
#endif
  array_free(&ls);
  schur->Q_ls = setup_Fxi_Q(&schur->A_ls, crs->lg + crs->is, crs->in, c, bfr);

  // Setup the preconditioner for the Schur complement matrix
  schur->M = schur_precond_setup(&schur->A_ll, &schur->A_ls, &schur->A_ss,
                                 &schur->A_sl, c, cr, bfr);

  return 0;
}

int schur_solve(scalar *x, scalar *b, struct coarse *crs, struct comm *c,
                buffer *bfr) {
  struct schur *schur = crs->solver;

  uint ln = crs->ln, in = crs->in, *idx = crs->idx;
  scalar *rhs = (scalar *)tcalloc(scalar, ln > in ? ln : in);
  scalar *zl = (scalar *)tcalloc(scalar, ln);
  scalar *xl = (scalar *)tcalloc(scalar, in + ln), *xi = xl + ln;

  // Solve: A_ll z_l = r_l
  for (uint i = 0; i < ln; i++)
    rhs[i] = b[idx[i]];
#ifdef DUMPI
  printf("%d bl = ", c->id);
  for (uint i = 0; i < ln; i++)
    printf("%lf ", rhs[i]);
  printf("\n");
  fflush(stdout);
  printf("%d bi = ", c->id);
  for (uint i = 0; i < in; i++)
    printf("%lf ", b[idx[ln + i]]);
  printf("\n");
  fflush(stdout);
#endif

  comm_barrier(c);
  double t = comm_time();
  cholesky_solve(zl, &schur->A_ll, rhs);
  t = comm_time() - t;
  if (c->id == 0)
    printf("\t cholesky_solve: %lf\n", t);

#ifdef DUMPI
  printf("%d zl = ", c->id);
  for (uint i = 0; i < ln; i++)
    printf("%lf ", zl[i]);
  printf("\n");
  fflush(stdout);
#endif

  // Solve: A_ss x_i = fi where fi = r_i - E zl
  Ezl(rhs, &schur->A_sl, schur->Q_sl, zl, crs->ls, in, bfr);
#ifdef DUMPI
  printf("%d Ezl = ", c->id);
  for (uint i = 0; i < in; i++)
    printf("%lf ", rhs[i]);
  printf("\n");
  fflush(stdout);
#endif

  for (uint i = 0; i < in; i++)
    rhs[i] = b[idx[ln + i]] - rhs[i];
#ifdef DUMPI
  printf("%d fi = ", c->id);
  for (uint i = 0; i < in; i++)
    printf("%lf ", rhs[i]);
  printf("\n");
  fflush(stdout);
#endif

  assert(in == schur->A_ss.rn);

  comm_barrier(c);
  t = comm_time();
  project(xi, rhs, schur, crs->ls, c, 100, 1, 1, bfr);
  t = comm_time() - t;
  if (c->id == 0)
    printf("\tproject: %lf\n", t);

#ifdef DUMPI
  printf("%d xi = ", c->id);
  for (uint i = 0; i < in; i++)
    printf("%lf ", xi[i]);
  printf("\n");
  fflush(stdout);
#endif

  // Solve A_ll xl = fl where fl = r_l - F xi
  for (uint i = 0; i < ln; i++)
    rhs[i] = 0;
  Fxi(rhs, &schur->A_ls, schur->Q_ls, xi, crs->ls, in, bfr);
#ifdef DUMPI
  printf("%d Fxi = ", c->id);
  for (uint i = 0; i < ln; i++)
    printf("%lf ", rhs[i]);
  printf("\n");
  fflush(stdout);
#endif

  for (uint i = 0; i < ln; i++)
    rhs[i] = b[idx[i]] - rhs[i];

  comm_barrier(c);
  t = comm_time();
  cholesky_solve(xl, &schur->A_ll, rhs);
  t = comm_time() - t;
  if (c->id == 0)
    printf("\tcholesky_solve 2: %lf\n", t);

#ifdef DUMPI
  printf("%d xl = ", c->id);
  for (uint i = 0; i < ln; i++)
    printf("%lf ", xl[i]);
  printf("\n");
  fflush(stdout);
#endif

  for (uint i = 0; i < ln + in; i++)
    x[idx[i]] = xl[i];

  free(rhs);
  free(zl);
  free(xl);

  return 0;
}

int schur_free(struct coarse *crs) {
  struct schur *schur = (struct schur *)crs->solver;
  if (schur != NULL) {
    mat_free(&schur->A_ll);
    par_mat_free(&schur->A_ls);
    if (schur->Q_ls != NULL)
      gs_free(schur->Q_ls), schur->Q_ls = NULL;
    par_mat_free(&schur->A_sl);
    if (schur->Q_sl != NULL)
      gs_free(schur->Q_sl), schur->Q_sl = NULL;
    par_mat_free(&schur->A_ss);
    if (schur->Q_ss != NULL)
      gs_free(schur->Q_ss), schur->Q_ss = NULL;
    if (schur->M != NULL)
      mg_free(schur->M), schur->M = NULL;
    free(schur), schur = NULL;
  }

  return 0;
}
