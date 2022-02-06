#include "coarse.h"
#include "multigrid.h"
#include <math.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define FREE(ptr, x)                                                           \
  {                                                                            \
    if (ptr->x != NULL)                                                        \
      free(ptr->x);                                                            \
  }

#define CSC 0
#define CSR 1

#define IS_CSC(A) ((A)->type == CSC)
#define IS_CSR(A) ((A)->type == CSR)
#define NO_DIAG(A) ((A)->diag_val == NULL)

//#define DUMPB
//#define DUMPE
//#define DUMPF
//#define DUMPS
//#define DUMPG
//#define DUMPW
//#define DUMPWG
//#define DUMPP

static void comm_split(const struct comm *old, const int bin, const int key,
                       struct comm *new_) {
  MPI_Comm new_comm;
  MPI_Comm_split(old->c, bin, key, &new_comm);
  comm_init(new_, new_comm);
  MPI_Comm_free(&new_comm);
}

//------------------------------------------------------------------------------
// Number rows, local first then interface. Returns global number of local
// elements.
static int number_rows(ulong *elem, ulong ng[2], const slong *vtx,
                       const uint nelt, const int nv, const struct comm *ci,
                       buffer *bfr) {
  uint npts = nelt * nv;
  int nnz = npts > 0;

  struct comm c;
  comm_split(ci, nnz, ci->id, &c);

  if (nnz == 0)
    return 1;

  int *dof = tcalloc(int, npts), level = 1;
  uint i, j;
  while (c.np > 1) {
    struct gs_data *gsh = gs_setup(vtx, npts, &c, 0, gs_pairwise, 0);

    int bin = c.id >= (c.np + 1) / 2;
    assert(bin == 0 || bin == 1);
    for (i = 0; i < npts; i++)
      dof[i] = bin;

    gs(dof, gs_int, gs_add, 0, gsh, bfr);

    if (bin == 1)
      for (i = 0; i < npts; i++)
        dof[i] = 0;

    gs(dof, gs_int, gs_add, 0, gsh, bfr);

    for (i = 0; i < nelt; i++) {
      for (j = 0; j < nv; j++) {
        if (dof[i * nv + j] > 0 && !elem[i]) {
          elem[i] = level;
          break;
        }
      }
    }

    gs_free(gsh);

    struct comm t;
    comm_split(&c, bin, c.id, &t);
    comm_free(&c);
    comm_dup(&c, &t);
    comm_free(&t);

    level++;
  }
  comm_free(&c);

  uint ni = 0;
  for (i = 0; i < nelt; i++)
    if (elem[i] > 0)
      ni++;
  uint nl = nelt - ni;

  slong in[2] = {nl, ni}, out[2][2], buf[2][2];
  comm_scan(out, ci, gs_long, gs_add, in, 2, buf);
  slong sl = out[0][0] + 1, si = out[0][1];
  slong gl = out[1][0] + 1;

  for (i = 0; i < nelt; i++)
    if (elem[i] > 0)
      elem[i] = gl + si++;
    else
      elem[i] = sl++;

  ng[0] = out[1][0], ng[1] = out[1][1];

  free(dof);
}

//------------------------------------------------------------------------------
// Find neighbors in the graph
//
struct nbr {
  ulong r, c;
  uint proc;
};

void find_nbrs(struct array *arr, const ulong *eid, const slong *vtx,
               const uint nelt, const int nv, const struct comm *c,
               struct crystal *cr, buffer *buf) {
  struct array vertices;
  array_init(struct nbr, &vertices, nelt * nv);

  struct nbr v;
  uint i, j;
  for (i = 0; i < nelt; i++) {
    v.r = eid[i];
    for (j = 0; j < nv; j++) {
      v.c = vtx[i * nv + j], v.proc = v.c % c->np;
      array_cat(struct nbr, &vertices, &v, 1);
    }
  }

  sarray_transfer(struct nbr, &vertices, proc, 1, cr);

  sarray_sort(struct nbr, vertices.ptr, vertices.n, c, 1, buf);
  struct nbr *vptr = (struct nbr *)vertices.ptr;
  uint vn = vertices.n;

  // FIXME: Assumes quads or hexes
  struct nbr t;
  uint s = 0, e;
  while (s < vn) {
    e = s + 1;
    while (e < vn && vptr[s].c == vptr[e].c)
      e++;
    for (i = s; i < e; i++) {
      t = vptr[i];
      for (j = s; j < e; j++) {
        t.c = vptr[j].r;
        array_cat(struct nbr, arr, &t, 1);
      }
    }
    s = e;
  }

  sarray_transfer(struct nbr, arr, proc, 1, cr);
  array_free(&vertices);
}

//------------------------------------------------------------------------------
// `mat` struct for local matrices
//
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

// Compress array by summing up entries which share the same (r, c) values
// and the array is modified in place. Also, the diagonal entries are modified
// to ensure all the
static int compress_nbrs(struct array *eij, struct array *nbr, buffer *bfr) {
  eij->n = 0;
  if (nbr->n == 0)
    return 1;

  sarray_sort_2(struct nbr, nbr->ptr, nbr->n, r, 1, c, 1, bfr);
  struct nbr *ptr = (struct nbr *)nbr->ptr;

  struct mat_ij m;
  m.idx = 0;

  sint i = 0;
  while (i < nbr->n) {
    m.r = ptr[i].r, m.c = ptr[i].c;

    sint j = i + 1;
    while (j < nbr->n && ptr[j].r == ptr[i].r && ptr[j].c == ptr[i].c)
      j++;

    m.v = i - j; // = - (j - i)
    array_cat(struct mat_ij, eij, &m, 1);
    i = j;
  }

  // Now make sure the row sum is zero
  struct mat_ij *pe = (struct mat_ij *)eij->ptr;
  i = 0;
  while (i < eij->n) {
    sint j = i, k = -1, s = 0;
    while (j < eij->n && pe[j].r == pe[i].r) {
      if (pe[j].r == pe[j].c)
        k = j;
      else
        s += pe[j].v;
      j++;
    }
    assert(k >= 0);
    pe[k].v = -s;
    i = j;
  }
}

static int csr_setup(struct mat *mat, struct array *entries, int sep,
                     buffer *buf) {
  uint nnz = entries->n;
  if (nnz == 0) {
    mat->start = mat->n = 0;
    mat->Lp = mat->Li = NULL;
    mat->L = mat->D = NULL;
    mat->col = NULL;
    return 0;
  }

  // Renumber cols and rows
  sarray_sort(struct mat_ij, entries->ptr, entries->n, c, 1, buf);
  struct mat_ij *ptr = (struct mat_ij *)entries->ptr;

  // Reserve enough memory for work arrays
  buffer_reserve(buf, sizeof(ulong) * nnz);
  ulong *cols = (ulong *)buf->ptr;

  ulong n = cols[0] = ptr[0].c;
  uint i;
  for (ptr[0].c = 0, i = 1; i < nnz; i++) {
    if (ptr[i].c == n)
      ptr[i].c = ptr[i - 1].c;
    else {
      n = ptr[i].c, ptr[i].c = ptr[i - 1].c + 1;
      cols[ptr[i].c] = n;
    }
  }
  uint nc = ptr[nnz - 1].c + 1;

  mat->col = tcalloc(ulong, nc);
  memcpy(mat->col, cols, sizeof(ulong) * nc);

  sarray_sort(struct mat_ij, entries->ptr, entries->n, r, 1, buf);
  ptr = (struct mat_ij *)entries->ptr;
  mat->start = ptr[0].r;

  for (n = ptr[0].r, ptr[0].r = 0, i = 1; i < nnz; i++)
    if (ptr[i].r == n)
      ptr[i].r = ptr[i - 1].r;
    else
      n = ptr[i].r, ptr[i].r = ptr[i - 1].r + 1;
  uint nr = ptr[nnz - 1].r + 1;

  // Reserve enough memory for work arrays
  buffer_reserve(buf, sizeof(struct mat_ij) * nnz);
  struct mat_ij *unique = (struct mat_ij *)buf->ptr;

  // Setup the matrix, separate diagonal
  mat->n = nr;
  uint *Lp = mat->Lp = (uint *)tcalloc(uint, nr + 1);

  sarray_sort_2(struct mat_ij, entries->ptr, entries->n, r, 1, c, 1, buf);
  ptr = (struct mat_ij *)entries->ptr;

  sep = sep > 0 ? 1 : 0;
  Lp[0] = 0, unique[0] = ptr[0];
  uint j;
  for (nr = 1, i = 1, j = 0; i < nnz; i++) {
    if ((unique[j].r != ptr[i].r) || (unique[j].c != ptr[i].c)) {
      if (unique[j].r != ptr[i].r)
        Lp[nr] = j + 1 - sep * nr, nr++;
      unique[++j] = ptr[i];
    } else
      unique[j].v += ptr[i].v;
  }
  Lp[mat->n] = ++j - sep * mat->n;

  mat->Li = (uint *)tcalloc(uint, j - sep * mat->n);
  mat->L = (scalar *)tcalloc(scalar, j - sep * mat->n);

  uint nadj;
  if (sep) {
    mat->D = (scalar *)tcalloc(scalar, mat->n);
    uint ndiag;
    for (i = ndiag = nadj = 0; i < j; i++) {
      if (mat->start + unique[i].r == mat->col[unique[i].c])
        mat->D[ndiag++] = unique[i].v;
      else
        mat->L[nadj] = unique[i].v, mat->Li[nadj++] = unique[i].c;
    }
  } else {
    mat->D = NULL;
    for (i = nadj = 0; i < j; i++)
      mat->L[nadj] = unique[i].v, mat->Li[nadj++] = unique[i].c;
  }

  return 0;
}

int mat_print(struct mat *mat) {
  uint i, j;
  for (i = 0; i < mat->n; i++) {
    for (j = mat->Lp[i]; j < mat->Lp[i + 1]; j++)
      printf("%ld %ld %lf\n", mat->start + i, mat->col[mat->Li[j]], mat->L[j]);
    if (mat->D != NULL)
      printf("%ld %ld %lf\n", mat->start + i, mat->start + i, mat->D[i]);
  }

  return 0;
}

static int mat_free(struct mat *mat) {
  FREE(mat, Lp);
  FREE(mat, Li);
  FREE(mat, L);
  FREE(mat, D);
  FREE(mat, col);
  return 0;
}

//------------------------------------------------------------------------------
// `par_mat` matrix for parallel distributed matrices
//
//
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
  uint *diag_idx;
  scalar *diag_val;
};

static int par_csr_setup(struct par_mat *mat, struct array *entries, int sd,
                         buffer *buf) {
  mat->type = CSR;
  if (entries == NULL || entries->n == 0) {
    mat->cn = mat->rn = 0;
    mat->diag_idx = mat->adj_off = mat->adj_idx = NULL;
    mat->adj_val = mat->diag_val = NULL;
    mat->rows = mat->cols = NULL;
    return 0;
  }

  sarray_sort(struct mat_ij, entries->ptr, entries->n, c, 1, buf);
  struct mat_ij *ptr = (struct mat_ij *)entries->ptr;

  // Reserve enough memory for work arrays
  uint nnz = entries->n;
  buffer_reserve(buf, (sizeof(struct mat_ij) + 2 * sizeof(ulong)) * (nnz + 1));

  ulong *cols = (ulong *)buf->ptr;
  cols[0] = ptr[0].c, ptr[0].idx = 0, mat->cn = 1;
  for (uint i = 1; i < nnz; i++) {
    if (ptr[i - 1].c != ptr[i].c)
      cols[mat->cn] = ptr[i].c, ptr[i].idx = mat->cn++;
    else
      ptr[i].idx = ptr[i - 1].idx;
  }

  mat->cols = (ulong *)tcalloc(ulong, mat->cn);
  memcpy(mat->cols, cols, sizeof(ulong) * mat->cn);

  sarray_sort_2(struct mat_ij, entries->ptr, entries->n, r, 1, c, 1, buf);
  ptr = entries->ptr;

  sd = sd != 0; // sd needs to be 1 or 0

  uint *adj_off = (uint *)buf->ptr;
  ulong *rows = (ulong *)(adj_off + nnz + 1);
  struct mat_ij *unique = (struct mat_ij *)(rows + (nnz + 1));

  adj_off[0] = 0, unique[0] = ptr[0], rows[0] = ptr[0].r;
  uint j = 0;
  for (uint i = mat->rn = 1; i < nnz; i++) {
    if ((unique[j].r != ptr[i].r) || (unique[j].c != ptr[i].c)) {
      if (unique[j].r != ptr[i].r) {
        adj_off[mat->rn] = j + 1 - sd * mat->rn;
        rows[mat->rn++] = ptr[i].r;
      }
      unique[++j] = ptr[i];
    } else
      unique[j].v += ptr[i].v;
  }
  adj_off[mat->rn] = ++j - sd * mat->rn;

  mat->rows = (ulong *)tcalloc(ulong, mat->rn);
  memcpy(mat->rows, rows, mat->rn * sizeof(ulong));

  mat->adj_off = (uint *)tcalloc(uint, mat->rn + 1);
  memcpy(mat->adj_off, adj_off, (mat->rn + 1) * sizeof(uint));
  mat->adj_idx = (uint *)tcalloc(uint, j - sd * mat->rn);
  mat->adj_val = (scalar *)tcalloc(scalar, j - sd * mat->rn);

  if (sd) {
    mat->diag_idx = (uint *)tcalloc(uint, mat->rn);
    mat->diag_val = (scalar *)tcalloc(scalar, mat->rn);
  } else {
    mat->diag_idx = NULL;
    mat->diag_val = NULL;
  }

  uint ndiag, nadj, i;
  for (i = ndiag = nadj = 0; i < j; i++) {
    if (unique[i].r == unique[i].c && sd) {
      mat->diag_idx[ndiag] = unique[i].idx;
      mat->diag_val[ndiag++] = unique[i].v;
    } else {
      mat->adj_idx[nadj] = unique[i].idx;
      mat->adj_val[nadj++] = unique[i].v;
    }
  }

  // Sanity check
  assert(ndiag == sd * mat->rn);
  assert(nadj == j - sd * mat->rn);

  return 0;
}

static int par_csc_setup(struct par_mat *mat, struct array *entries, int sd,
                         buffer *buf) {
  mat->type = CSC;
  if (entries == NULL || entries->n == 0) {
    mat->cn = mat->rn = 0;
    mat->diag_idx = mat->adj_off = mat->adj_idx = NULL;
    mat->adj_val = mat->diag_val = NULL;
    mat->rows = mat->cols = NULL;
    return 0;
  }

  sarray_sort(struct mat_ij, entries->ptr, entries->n, r, 1, buf);
  struct mat_ij *ptr = (struct mat_ij *)entries->ptr;

  // Reserve enough memory for work arrays
  uint nnz = entries->n;
  buffer_reserve(buf, (sizeof(struct mat_ij) + 2 * sizeof(ulong)) * (nnz + 1));

  ulong *rows = (ulong *)buf->ptr;
  rows[0] = ptr[0].r, ptr[0].idx = 0, mat->rn = 1;
  for (uint i = 1; i < nnz; i++) {
    if (ptr[i - 1].r != ptr[i].r)
      rows[mat->rn] = ptr[i].r, ptr[i].idx = mat->rn++;
    else
      ptr[i].idx = ptr[i - 1].idx;
  }

  mat->rows = tcalloc(ulong, mat->rn);
  memcpy(mat->rows, rows, sizeof(ulong) * mat->rn);

  sarray_sort_2(struct mat_ij, entries->ptr, entries->n, c, 1, r, 1, buf);
  ptr = entries->ptr;

  sd = sd != 0; // sd needs to be 1 or 0

  uint *adj_off = (uint *)buf->ptr;
  ulong *cols = (ulong *)(adj_off + nnz + 1);
  struct mat_ij *unique = (struct mat_ij *)(cols + (nnz + 1));

  adj_off[0] = 0, unique[0] = ptr[0], cols[0] = ptr[0].c;
  uint j = 0;
  for (uint i = mat->cn = 1; i < nnz; i++) {
    if ((unique[j].r != ptr[i].r) || (unique[j].c != ptr[i].c)) {
      if (unique[j].c != ptr[i].c) {
        adj_off[mat->cn] = j + 1 - sd * mat->cn;
        cols[mat->cn++] = ptr[i].c;
      }
      unique[++j] = ptr[i];
    } else
      unique[j].v += ptr[i].v;
  }
  adj_off[mat->cn] = ++j - sd * mat->cn;

  mat->cols = (ulong *)tcalloc(ulong, mat->cn);
  memcpy(mat->cols, cols, mat->cn * sizeof(ulong));

  mat->adj_off = (uint *)tcalloc(uint, mat->cn + 1);
  memcpy(mat->adj_off, adj_off, (mat->cn + 1) * sizeof(uint));
  mat->adj_idx = (uint *)tcalloc(uint, j - sd * mat->cn);
  mat->adj_val = (scalar *)tcalloc(scalar, j - sd * mat->cn);

  if (sd) {
    mat->diag_idx = (uint *)tcalloc(uint, mat->cn);
    mat->diag_val = (scalar *)tcalloc(scalar, mat->cn);
  } else {
    mat->diag_idx = NULL;
    mat->diag_val = NULL;
  }

  uint ndiag, nadj, i;
  for (i = ndiag = nadj = 0; i < j; i++) {
    if (unique[i].r == unique[i].c && sd) {
      mat->diag_idx[ndiag] = unique[i].idx;
      mat->diag_val[ndiag++] = unique[i].v;
    } else {
      mat->adj_idx[nadj] = unique[i].idx;
      mat->adj_val[nadj++] = unique[i].v;
    }
  }

  // Sanity check
  assert(ndiag == sd * mat->cn);
  assert(nadj == j - sd * mat->cn);

  return 0;
}

static void par_mat_print(struct par_mat *A) {
  uint i, j;
  if (IS_CSR(A)) {
    for (i = 0; i < A->rn; i++) {
      for (j = A->adj_off[i]; j < A->adj_off[i + 1]; j++)
        printf("%ld %ld %lf\n", A->rows[i], A->cols[A->adj_idx[j]],
               A->adj_val[j]);
      if (A->diag_val != NULL)
        printf("%ld %ld %lf\n", A->rows[i], A->rows[i], A->diag_val[i]);
    }
  } else if (IS_CSC(A)) {
    for (i = 0; i < A->cn; i++) {
      for (j = A->adj_off[i]; j < A->adj_off[i + 1]; j++)
        printf("%ld %ld %lf\n", A->rows[A->adj_idx[j]], A->cols[i],
               A->adj_val[j]);
      if (A->diag_val != NULL)
        printf("%ld %ld %lf\n", A->cols[i], A->cols[i], A->diag_val[i]);
    }
  }
}

int par_mat_free(struct par_mat *A) {
  FREE(A, rows);
  FREE(A, cols);
  FREE(A, adj_off);
  FREE(A, adj_idx);
  FREE(A, adj_val);
  FREE(A, diag_idx);
  FREE(A, diag_val);

  return 0;
}

struct par_mat *par_csr_setup_con(const uint nelt, const ulong *eid,
                                  const slong *vtx, int nv, int sep,
                                  struct comm *c, struct crystal *cr,
                                  buffer *bfr) {
  struct array nbrs, eij;
  array_init(struct nbr, &nbrs, 100);
  array_init(struct mat_ij, &eij, 100);

  find_nbrs(&nbrs, eid, vtx, nelt, nv, c, cr, bfr);
  compress_nbrs(&eij, &nbrs, bfr);
  array_free(&nbrs);

  struct par_mat *M = tcalloc(struct par_mat, 1);
  par_csr_setup(M, &eij, sep, bfr);
  array_free(&eij);

  return M;
}

static int compress_mat_ij(struct array *eij, struct array *entries,
                           buffer *bfr) {
  eij->n = 0;
  if (entries->n == 0)
    return 1;

  sarray_sort_2(struct mat_ij, entries->ptr, entries->n, r, 1, c, 1, bfr);
  struct mat_ij *ptr = (struct mat_ij *)entries->ptr;

  struct mat_ij m;
  m.idx = 0;

  uint i = 0;
  while (i < entries->n) {
    m = ptr[i];
    uint j = i + 1;
    while (j < entries->n && ptr[j].r == ptr[i].r && ptr[j].c == ptr[i].c)
      m.v += ptr[j].v, j++;

    array_cat(struct mat_ij, eij, &m, 1);
    i = j;
  }

  // Now make sure the row sum is zero
  struct mat_ij *pe = (struct mat_ij *)eij->ptr;
  i = 0;
  while (i < eij->n) {
    sint j = i, k = -1, s = 0;
    while (j < eij->n && pe[j].r == pe[i].r) {
#if 0
      if (eij->n < 5)
        printf("r = %lld c = %lld\n", pe[j].r, pe[j].c);
#endif
      if (pe[j].r == pe[j].c)
        k = j;
      else
        s += pe[j].v;
      j++;
    }
#if 0
    printf("ii = %d n = %d\n", i, eij->n);
#endif
    assert(k >= 0);
    pe[k].v = -s;
    i = j;
  }
}

struct par_mat *par_csr_setup_ext(struct array *entries, int sep, buffer *bfr) {
  struct array eij;
  array_init(struct mat_ij, &eij, 100);

  compress_mat_ij(&eij, entries, bfr);

  struct par_mat *M = tcalloc(struct par_mat, 1);
  par_csr_setup(M, &eij, sep, bfr);

  array_free(&eij);

  return M;
}

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

//-----------------------------------------------------------------------------
// Schur setup
//
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
static int schur_setup_G(struct par_mat *G, const struct mat *L,
                         const struct par_mat *F, const ulong *srows,
                         const uint srn, const struct comm *c,
                         struct crystal *const cr, buffer *bfr) {
  assert(IS_CSR(F));
  assert(NO_DIAG(F));

  uint n = L->n;
  scalar *b = tcalloc(scalar, 2 * n);
  scalar *x = b + n;

  struct array gij;
  array_init(struct mat_ij, &gij, 100);

  struct mat_ij m;
  uint i, j, k, ke;
  for (i = 0; i < F->rn; i++) {
    b[F->rows[i] - L->start] = 1;
    cholesky_lower_solve(x, L, b);
#if 0
    printf("L = ");
    for (j = 0; j < n; j++)
      printf("%lf ", x[j]);
    printf("\n");
    fflush(stdout);
#endif

    // Calculate F: i^th row of F is multiplied by each element of i^th
    // column of L_B^-1
    for (k = F->adj_off[i], ke = F->adj_off[i + 1]; k < ke; k++) {
      m.c = F->cols[F->adj_idx[k]];
      for (j = 0; j < n; j++) {
        m.r = L->start + j;
        m.v = F->adj_val[k] * x[j];
        array_cat(struct mat_ij, &gij, &m, 1);
      }
    }

    b[F->rows[i] - L->start] = 0;
    for (j = 0; j < n; j++)
      x[j] = 0;
  }

  sarray_sort_2(struct mat_ij, gij.ptr, gij.n, r, 1, c, 1, bfr);
  struct mat_ij *ptr = (struct mat_ij *)gij.ptr;

  buffer_reserve(bfr, (sizeof(slong) + sizeof(sint)) * (gij.n + 1));
  slong *cols = (slong *)bfr->ptr;
  sint *owners = (sint *)(cols + gij.n + 1);

  struct array unique;
  array_init(struct mat_ij, &unique, 100);
  for (i = k = 0; i < gij.n; k++) {
    for (j = i + 1; j < gij.n && ptr[j].r == ptr[i].r && ptr[j].c == ptr[i].c;
         j++)
      ptr[i].v += ptr[j].v;

    array_cat(struct mat_ij, &unique, &ptr[i], 1);
    cols[k] = ptr[i].c;
    owners[k] = S_owns_row(cols[k], srows, srn) < srn ? c->id : -1;
    i = j;
  }
  array_free(&gij);
  assert(k == unique.n); // Sanity check

  // Set the destination processor for elements in unique. Since W will be in
  // CSR and share the same row distribution as S, we use gslib and row
  // distribution of S to find the destination processor.
  struct gs_data *gsh = gs_setup(cols, k, c, 0, gs_pairwise, 1);
  gs(owners, gs_int, gs_max, 0, gsh, bfr);
  free(gsh);

  ptr = unique.ptr;
  for (i = 0; i < unique.n; i++)
    ptr[i].p = owners[i];

  sarray_transfer(struct mat_ij, &unique, p, 0, cr);

  par_csc_setup(G, &unique, 0, bfr);
#ifdef DUMPG
  par_mat_print(G);
#endif

  array_free(&unique);
  free(b);

  return 0;
}

// Calculate W = E x U_{B}^{-1} where B = L_{B} U_{B}. E is in CSC format.
// W will be in CSR format and distributed by rows similar to distribution of S.
static int schur_setup_W(struct par_mat *W, const struct mat *L,
                         const struct par_mat *Er, const ulong *srows,
                         const uint srn, const struct comm *c,
                         struct crystal *const cr, buffer *bfr) {
  assert(IS_CSR(Er));
  assert(NO_DIAG(Er));

  uint nnz = 0;
  if (Er->rn > 0)
    nnz = Er->adj_off[Er->rn];

  // Setup E as a CSC matrix. First, find the owner processor of each column
  // with gslib. Then send the matrix entries to relevant processor and setup
  // E as a CSC matrix.
  buffer_reserve(bfr, (sizeof(sint) + sizeof(slong)) * nnz);
  slong *cols = (slong *)bfr->ptr;
  sint *owners = (sint *)(cols + nnz);

  uint i, k, ke;
  for (i = 0; i < Er->rn; i++) {
    for (k = Er->adj_off[i], ke = Er->adj_off[i + 1]; k != ke; k++) {
      cols[k] = Er->cols[Er->adj_idx[k]];
      if (cols[k] >= L->start && cols[k] < L->start + L->n)
        owners[k] = c->id;
      else
        owners[k] = -1;
    }
  }

  struct gs_data *gsh = gs_setup(cols, nnz, c, 0, gs_pairwise, 0);
  gs(owners, gs_int, gs_max, 0, gsh, bfr);
  gs_free(gsh);

  struct array eij;
  array_init(struct mat_ij, &eij, 100);
  struct mat_ij *ptr = (struct mat_ij *)eij.ptr;

  struct mat_ij m;
  for (i = 0; i < Er->rn; i++) {
    m.r = Er->rows[i];
    for (k = Er->adj_off[i], ke = Er->adj_off[i + 1]; k != ke; k++) {
      m.c = cols[k], m.p = owners[k], m.v = Er->adj_val[k];
      array_cat(struct mat_ij, &eij, &m, 1);
    }
  }

  sarray_transfer(struct mat_ij, &eij, p, 1, cr);

  struct par_mat E;
  par_csc_setup(&E, &eij, 0, bfr);
  array_free(&eij);

  // Multiply E by U_B^{-1} now. Columns of U_B^{-1} are found one by one and
  // then E is multiplied by each column.
  uint n = L->n;
  scalar *b = tcalloc(scalar, 2 * n);
  scalar *x = b + n;

  struct array wij;
  array_init(struct mat_ij, &wij, 100);

  uint j;
  for (i = 0; i < n; i++) {
    b[i] = 1;
    cholesky_upper_solve(x, L, b);
#if 0
    printf("U = ");
    for (j = 0; j < n; j++)
      printf("%lf ", x[j]);
    printf("\n");
    fflush(stdout);
#endif

    // Multiply E by x: i^th col of E is multiplied by element x[i]
    for (j = 0; j < E.cn; j++) {
      m.c = L->start + i;
      for (k = E.adj_off[j], ke = E.adj_off[j + 1]; k < ke; k++) {
        m.r = E.rows[E.adj_idx[k]];
        m.v = E.adj_val[k] * x[E.cols[j] - L->start];
        array_cat(struct mat_ij, &wij, &m, 1);
      }
    }

    b[i] = 0;
    for (j = 0; j < n; j++)
      x[j] = 0;
  }
  par_mat_free(&E);

  sarray_sort_2(struct mat_ij, wij.ptr, wij.n, r, 1, c, 1, bfr);
  ptr = (struct mat_ij *)wij.ptr;

  buffer_reserve(bfr, (sizeof(sint) + sizeof(slong)) * (wij.n + 1));
  slong *rows = (slong *)bfr->ptr;
  owners = (sint *)(rows + wij.n + 1);

  struct array unique;
  array_init(struct mat_ij, &unique, 100);
  for (i = k = 0; i < wij.n; k++) {
    for (j = i + 1; j < wij.n && ptr[j].r == ptr[i].r && ptr[j].c == ptr[i].c;
         j++)
      ptr[i].v += ptr[j].v;

    array_cat(struct mat_ij, &unique, &ptr[i], 1);
    rows[k] = ptr[i].r;
    owners[k] = S_owns_row(rows[k], srows, srn) < srn ? c->id : -1;
    i = j;
  }
  array_free(&wij);
  assert(k == unique.n); // Sanity check

  // Set the destination processor for elements in unique. Since W will be in
  // CSR and share the same row distribution as S, we use gslib and row
  // distribution of S to find the destination processor.
  gsh = gs_setup(rows, k, c, 0, gs_pairwise, 0);
  gs(owners, gs_int, gs_max, 0, gsh, bfr);
  free(gsh);

  ptr = unique.ptr;
  for (i = 0; i < unique.n; i++)
    ptr[i].p = owners[i];

  sarray_transfer(struct mat_ij, &unique, p, 0, cr);

  par_csr_setup(W, &unique, 0, bfr);
#ifdef DUMPW
  par_mat_print(W);
#endif

  array_free(&unique);
  free(b);

  return 0;
}

// C = A - B; A and B should be in CSR format with the same row
// distribution
static int sparse_sub(struct par_mat *C, const struct par_mat *A,
                      const struct par_mat *B, const struct comm *c,
                      struct crystal *const cr, buffer *bfr) {
  assert(IS_CSR(A));
  assert(IS_CSR(B));

  struct array cij;
  array_init(struct mat_ij, &cij, 100);

  struct mat_ij m;
  uint r, j, je;
  for (r = 0; r < A->rn; r++) {
    m.r = A->rows[r];
    for (j = A->adj_off[r], je = A->adj_off[r + 1]; j != je; j++) {
      m.c = A->cols[A->adj_idx[j]], m.v = A->adj_val[j];
      array_cat(struct mat_ij, &cij, &m, 1);
    }
  }

  for (r = 0; r < B->rn; r++) {
    m.r = B->rows[r];
    for (j = B->adj_off[r], je = B->adj_off[r + 1]; j != je; j++) {
      m.c = B->cols[B->adj_idx[j]], m.v = -B->adj_val[j];
      array_cat(struct mat_ij, &cij, &m, 1);
    }
  }

  sarray_sort_2(struct mat_ij, cij.ptr, cij.n, r, 1, c, 1, bfr);
  struct mat_ij *ptr = (struct mat_ij *)cij.ptr;

  struct array unique;
  array_init(struct mat_ij, &unique, 100);

  if (cij.n > 0) {
    uint i = 0;
    while (i < cij.n) {
      j = i;
      scalar s = 0;
      for (; j < cij.n && ptr[j].r == ptr[i].r && ptr[j].c == ptr[i].c; j++)
        s += ptr[j].v;
      m = ptr[i], m.v = s;
      array_cat(struct mat_ij, &unique, &m, 1);
      i = j;
    }
  }
  array_free(&cij);

  par_csr_setup(C, &unique, 1, bfr);
  array_free(&unique);

  return 0;
}

static int sparse_gemm(struct par_mat *S, const struct par_mat *W,
                       const struct par_mat *G, struct crystal *cr,
                       const struct comm *c, buffer *bfr) {
  // W is in CSR, G is in CSC; we multiply rows of W by shifting
  // the columns of G from processor to processor. This is not scalable
  // at all -- need to do a 2D partition of the matrices W and G.
  assert(IS_CSR(W));
  assert(IS_CSC(G));

  // Put G into an array to transfer from processor to processor
  struct array gij, sij;
  array_init(struct mat_ij, &gij, 100);
  array_init(struct mat_ij, &sij, 100);

  struct mat_ij m;
  uint i, j, je;
  for (i = 0; i < G->cn; i++) {
    m.c = G->cols[i], m.p = c->id;
    for (j = G->adj_off[i], je = G->adj_off[i + 1]; j != je; j++) {
      m.r = G->rows[G->adj_idx[j]];
      m.v = G->adj_val[j];
      array_cat(struct mat_ij, &gij, &m, 1);
    }
  }

  sarray_sort_2(struct mat_ij, gij.ptr, gij.n, r, 1, c, 1, bfr);
  struct mat_ij *ptr = (struct mat_ij *)gij.ptr;
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
          array_cat(struct mat_ij, &sij, &m, 1);
        }
      }
    }

    sint next = (c->id + 1) % c->np;
    for (i = 0; i < gij.n; i++)
      ptr[i].p = next;
    sarray_transfer(struct mat_ij, &gij, p, 0, cr);

    sarray_sort(struct mat_ij, gij.ptr, gij.n, idx, 0, bfr);
    ptr = gij.ptr;
  }

  par_csr_setup(S, &sij, 0, bfr);

  array_free(&gij);
  array_free(&sij);

  return 0;
}

static struct mg_data *
schur_precond_setup(const struct mat *L, const struct par_mat *F,
                    const struct par_mat *S, const struct par_mat *E,
                    const struct comm *c, struct crystal *cr, buffer *bfr) {
  // TODO: Sparsify W and G when they are built
  struct par_mat W, G, WG;
  schur_setup_G(&G, L, F, S->rows, S->rn, c, cr, bfr);
  schur_setup_W(&W, L, E, S->rows, S->rn, c, cr, bfr);
  sparse_gemm(&WG, &W, &G, cr, c, bfr);
#ifdef DUMPWG
  par_mat_print(&WG);
#endif

  // P is CSR
  struct par_mat *P = tcalloc(struct par_mat, 1);
  sparse_sub(P, S, &WG, c, cr, bfr);
#ifdef DUMPP
  par_mat_print(&P);
#endif

  par_mat_free(&W);
  par_mat_free(&G);
  par_mat_free(&WG);

  return mg_setup(P, 2, c, cr, bfr);
}

static struct mg_data *schur_setup(struct mat *A_ll, struct par_mat *A_ls,
                                   struct par_mat *A_sl, struct par_mat *A_ss,
                                   struct array *eij, const ulong ng[2],
                                   struct comm *c, struct crystal *cr,
                                   buffer *bfr) {
  // Setup A_ll
  struct array ll, ls, sl, ss;
  array_init(struct mat_ij, &ll, eij->n / 4 + 1);
  array_init(struct mat_ij, &ls, eij->n / 4 + 1);
  array_init(struct mat_ij, &sl, eij->n / 4 + 1);
  array_init(struct mat_ij, &ss, eij->n / 4 + 1);

  struct mat_ij *ptr = (struct mat_ij *)eij->ptr;
  for (uint i = 0; i < eij->n; i++) {
    if (ptr[i].r < ng[0])
      if (ptr[i].c < ng[0])
        array_cat(struct mat_ij, &ll, &ptr[i], 1);
      else
        array_cat(struct mat_ij, &sl, &ptr[i], 1);
    else if (ptr[i].c < ng[0])
      array_cat(struct mat_ij, &ls, &ptr[i], 1);
    else
      array_cat(struct mat_ij, &ss, &ptr[i], 1);
  }

  // Setup local block diagonal (B)
  struct mat All;
  csr_setup(&All, &ll, 0, bfr);
  cholesky_factor(A_ll, &All, bfr);
#ifdef DUMPB
  mat_print(&A_ll);
#endif
  mat_free(&All);

  // Setup E
  par_csr_setup(A_ls, &ls, 0, bfr);
#ifdef DUMPE
  par_mat_print(A_ls);
#endif

  // Setup F
  par_csr_setup(A_sl, &sl, 0, bfr);
#ifdef DUMPF
  par_mat_print(A_sl);
#endif

  // Setup S
  par_csr_setup(A_ss, &ss, 0, bfr);
#ifdef DUMPS
  par_mat_print(A_ss);
#endif

  array_free(&ll);
  array_free(&ls);
  array_free(&sl);
  array_free(&ss);

  // Setup the preconditioner for the Schur complement matrix
  return schur_precond_setup(A_ll, A_sl, A_ss, A_ls, c, cr, bfr);
}

//------------------------------------------------------------------------------
// Setup coarse grid system
// A_ll: Local dof of a processor (block diagonal across processors)
// A_sl (= A_ls^T): shared - local matrix
// A_ss: Shared dof freedom (matrix is split row wise)
//     |A_ll (B)  A_ls (F)|
//  A= |                  |
//     |A_sl (E)  A_ss (S)|
struct coarse {
  struct par_mat A_ls, A_ss, A_sl;
  struct mat A_ll;
  struct mg_data *M;
};

struct coarse *coarse_setup(uint nelt, int nv, slong const *vtx, struct comm *c,
                            buffer *bfr) {
  struct coarse *crs = tcalloc(struct coarse, 1);

  ulong *eid = tcalloc(ulong, nelt), ng[2];
  number_rows(eid, ng, vtx, nelt, nv, c, bfr);

  struct crystal cr;
  crystal_init(&cr, c);

  struct array nbrs;
  array_init(struct nbr, &nbrs, nelt);
  find_nbrs(&nbrs, eid, vtx, nelt, nv, c, &cr, bfr);
  free(eid);

  // convert `struct nbr` -> `struct mat_ij` and compress
  // entries which share the same (r, c) values. Set the
  // diagonal element to have zero row sum

  // Be a bit smarter than just using 100?
  struct array eij;
  array_init(struct mat_ij, &eij, 100);
  compress_nbrs(&eij, &nbrs, bfr);
  array_free(&nbrs);

  crs->M = schur_setup(&crs->A_ll, &crs->A_ls, &crs->A_sl, &crs->A_ss, &eij, ng,
                       c, &cr, bfr);
  array_free(&eij);
  crystal_free(&cr);

  return crs;
}

int coarse_solve(scalar *x, struct coarse *crs, scalar *b) {
  cholesky_solve(x, &crs->A_ll, b);
  return 0;
}

int coarse_free(struct coarse *crs) {
  mat_free(&crs->A_ll);
  par_mat_free(&crs->A_ls);
  par_mat_free(&crs->A_sl);
  par_mat_free(&crs->A_ss);
  mg_free(crs->M);
  free(crs);
  return 0;
}

#undef NO_DIAG
#undef IS_CSC
#undef IS_CSR
#undef CSC
#undef CSR
#undef MIN
#undef FREE
