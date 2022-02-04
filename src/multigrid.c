#include "multigrid.h"
#include "coarse.h"
#include <math.h>

struct mg_lvl {
  int nsmooth;
  scalar sigma;
  struct gs_data *J; // Interpolation from level i to i + 1
  struct gs_data *Q; // gs handle for matrix vector product
  struct par_mat *M;
};

struct mg_data {
  int nlevels;
  struct mg_lvl **levels;
  uint *level_off;
  scalar *buf;
};

// Following struct definitions should match what is defined in coarse.c
struct mat_ij {
  ulong r, c;
  uint idx, p;
  scalar v;
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
  uint *diag_idx;
  scalar *diag_val;
};

static int logbll(slong n, int a) {
  assert(a > 0);

  int k = 0;
  while (n > 1)
    n /= a, k++;

  return k;
}

static struct gs_data *setup_Q(const struct par_mat *M, const struct comm *c,
                               buffer *bfr) {
  // Setup gs handle for the mat-vec
  uint nnz = M->rn > 0 ? M->adj_off[M->rn] + (M->diag_val != NULL) * M->rn : 0;
  buffer_reserve(bfr, nnz * sizeof(slong));

  slong *ids = (slong *)bfr->ptr;
  uint i, j;
  for (i = 0; i < M->rn; i++)
    for (j = M->adj_off[i]; j < M->adj_off[i + 1]; j++)
      ids[j] = -M->cols[M->adj_idx[j]];
  for (i = 0; i < M->rn; i++)
    ids[j++] = M->rows[i];

  return gs_setup(ids, nnz, c, 0, gs_crystal_router, 0);
}

static void mg_lvl_setup(struct mg_data *d, const uint lvl, const int factor,
                         const struct comm *c, struct crystal *cr,
                         struct array *entries, buffer *bfr) {
  assert(lvl > 0);
  struct par_mat *M = d->levels[lvl - 1]->M;
  uint nnz = M->rn > 0 ? M->adj_off[M->rn] + (M->diag_val != NULL) * M->rn : 0;

  // Reserve enough space in the array
  array_reserve(struct mat_ij, entries, nnz);
  entries->n = 0;

  // Reserve enough memory for gs ids for interpolation
  slong *ids = (slong *)tcalloc(slong, 2 * M->rn);

  // Calculate coarse level parameters: ngc, npc, nelt, nrem
  slong out[2][1], bf[2][1], in = M->rn;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bf);
  ulong ng = out[1][0], ngc = (ng + factor - 1) / factor;
  uint npc = ngc < c->np ? ngc : c->np;
  uint nelt = ngc / npc, nrem = ngc % npc;

#define SETP(p, r, np, nelt, nrem)                                             \
  {                                                                            \
    if (nrem == 0)                                                             \
      p = (r - 1) / nelt;                                                      \
    else {                                                                     \
      uint s = np - nrem;                                                      \
      ulong t = nelt * s;                                                      \
      if (r <= t)                                                              \
        p = (r - 1) / nelt;                                                    \
      else                                                                     \
        p = s + (r - t) / (nelt + 1);                                          \
    }                                                                          \
  }

  uint i, j, k = 0, n = 0;
  struct mat_ij m;
  for (i = 0; i < M->rn; i++) {
    ulong r = M->rows[i];
    for (j = M->adj_off[i]; j < M->adj_off[i + 1]; j++) {
      m.r = (r + factor - 1) / factor;
      if (m.r > ngc)
        m.r--;
      m.c = (M->cols[M->adj_idx[j]] + factor - 1) / factor;
      if (m.c > ngc)
        m.c--;
      m.v = M->adj_val[j];
      assert(m.r > 0);
      assert(m.c > 0);
      SETP(m.p, m.r, npc, nelt, nrem);
      array_cat(struct mat_ij, entries, &m, 1);
    }
    m.r = (r + factor - 1) / factor;
    if (m.r > ngc)
      m.r--;
    m.c = m.r;
    assert(m.r > 0);
    assert(m.c > 0);
    m.v = M->diag_val[i];
    SETP(m.p, m.r, npc, nelt, nrem);
    ids[k++] = -m.c;
    array_cat(struct mat_ij, entries, &m, 1);
  }

#undef SETP

  sarray_transfer(struct mat_ij, entries, p, 0, cr);
  sarray_sort_2(struct mat_ij, entries->ptr, entries->n, r, 1, c, 1, bfr);

  struct mg_lvl *l = d->levels[lvl] = tcalloc(struct mg_lvl, 1);
  M = l->M = par_csr_setup_ext(entries, 1, bfr);

  // Setup gs ids for coarse level (rhs interpolation )
  ids = trealloc(slong, ids, k + M->rn);
  for (i = 0; i < M->rn; i++)
    ids[k++] = M->rows[i];

  d->levels[lvl - 1]->J = gs_setup(ids, k, c, 0, gs_pairwise, 0);
  free(ids);
  d->levels[lvl]->Q = setup_Q(M, c, bfr);
}

struct mg_data *mg_setup(const struct par_mat *M, const int factor,
                         const struct comm *c, struct crystal *cr,
                         buffer *bfr) {
  slong out[2][1], bf[2][1], in = M->rn;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bf);
  slong nelg = out[1][0];

  struct mg_data *d = tcalloc(struct mg_data, 1);
  d->nlevels = logbll(nelg, factor) + 1;
  d->levels = tcalloc(struct mg_lvl *, d->nlevels);
#if 0
  if (c->id == 0)
    printf("nelg = %lld factor = %d nlevels = %d\n", nelg, factor, d->nlevels);
#endif

  d->levels[0] = tcalloc(struct mg_lvl, 1);
  d->levels[0]->nsmooth = 2;
  d->levels[0]->sigma = 0.6;
  d->levels[0]->M = (struct par_mat *)M;
  d->levels[0]->Q = setup_Q(M, c, bfr);

  d->level_off = tcalloc(uint, d->nlevels + 1);
  d->level_off[0] = 0;
  d->level_off[1] = M->rn;

  uint nnz = M->rn > 0 ? M->adj_off[M->rn] + M->rn : 0;
  struct array entries;
  array_init(struct mat_ij, &entries, nnz);

  uint i;
  for (i = 1; i < d->nlevels; i++) {
    mg_lvl_setup(d, i, factor, c, cr, &entries, bfr);
    struct par_mat *M = d->levels[i]->M;
#if 0
    if (c->id == 0)
      printf("i = %d\n", i);
    mat_print(d->levels[i]->M);
    fflush(stdout);
#endif
    if (M->rn > 0 && M->adj_off[M->rn] + M->rn > nnz)
      nnz = M->adj_off[M->rn] + M->rn;
    d->level_off[i + 1] = d->level_off[i] + M->rn;
    d->levels[i]->nsmooth = 2;
    d->levels[i]->sigma = 0.6;
  }

  d->levels[i - 1]->J = NULL;
  d->buf = tcalloc(scalar, 5 * d->level_off[d->nlevels] + nnz);

  array_free(&entries);

  return d;
}

static void mat_vec(scalar *y, struct par_mat *M, scalar *x,
                    struct gs_data *gsh, scalar *buf, buffer *bfr) {
  uint n = M->rn > 0 ? M->rn : 0;
  uint nnz = n > 0 ? M->adj_off[n] : 0;
  uint i;
  for (i = 0; i < nnz; i++)
    buf[i] = 0.0; // Is this really necessary?
  for (i = 0; i < M->rn; i++)
    y[i] = buf[nnz + i] = x[i];
  gs(buf, gs_double, gs_add, 0, gsh, bfr);

  uint j, je;
  for (i = 0; i < n; i++) {
    for (y[i] *= M->diag_val[i], j = M->adj_off[i], je = M->adj_off[i + 1];
         j != je; j++)
      y[i] += M->adj_val[j] * buf[j];
  }
}

void mg_vcycle(scalar *u1, scalar *rhs, struct mg_data *d, struct comm *c,
               buffer *bfr) {
  int nlevels = d->nlevels;
  uint *lvl_off = d->level_off, nnz = lvl_off[nlevels];
  scalar *s = d->buf, *Gs = s + nnz, *r = Gs + nnz, *u = r + nnz;
  uint i;
  for (i = 0; i < nnz; i++)
    s[i] = Gs[i] = r[i] = u[i] = 0.0;
  for (i = 0; i < lvl_off[1]; i++)
    r[i] = rhs[i];

  scalar *buf = u + nnz;
  uint n, j, off;
  struct mg_lvl **lvls = d->levels;
  for (int lvl = 0; lvl < nlevels - 1; lvl++) {
    off = lvl_off[lvl], n = lvl_off[lvl + 1] - off;
    struct mg_lvl *l = lvls[lvl];
    struct par_mat *M = l->M;

    // u = sigma*inv(D)*rhs
    scalar sigma = l->sigma;
    for (j = 0; j < n; j++)
      u[off + j] = sigma * r[off + j] / M->diag_val[j];

    // G*u
    mat_vec(Gs + off, M, u + off, l->Q, buf, bfr);

    // r = rhs - Gu
    for (j = 0; j < n; j++)
      r[off + j] = r[off + j] - Gs[off + j];

    for (i = 0; i < l->nsmooth; i++) {
      sigma = sigma + 0.066666 / l->nsmooth;

      // s = sigma*inv(D)*r, u = u + s
      for (j = 0; j < n; j++) {
        s[off + j] = sigma * r[off + j] / M->diag_val[j];
        u[off + j] += s[off + j];
      }

      // G*s
      mat_vec(Gs + off, M, s + off, l->Q, buf, bfr);

      // r=r-Gs
      for (j = 0; j < n; j++)
        r[off + j] = r[off + j] - Gs[off + j];
    }

    // Interpolate to coarser level
    gs(r + off, gs_double, gs_add, 1, l->J, bfr);
  }

  // Coarsest level
  off = lvl_off[nlevels - 1], n = lvl_off[nlevels] - off;

  if (n == 1) {
    struct mg_lvl *l = lvls[nlevels - 1];
    struct par_mat *M = l->M;
    assert(M->rn == 1);
    if (fabs(M->diag_val[0]) > sqrt(TOL))
      u[off] = r[off] / M->diag_val[0];
    else
      u[off] = 0.0;
    r[off] = u[off];
  }

  scalar over = 1.33333;
  for (int lvl = nlevels - 2; lvl >= 0; lvl--) {
    struct mg_lvl *l = lvls[lvl];
    off = lvl_off[lvl];
    // J*e
    gs(r + off, gs_double, gs_add, 0, l->J, bfr);

    // u = u + over*J*e
    n = lvl_off[lvl + 1] - off;
    for (j = 0; j < n; j++)
      r[off + j] = over * r[off + j] + u[off + j];
  }

  // Avoid this
  for (i = 0; i < lvl_off[1]; i++)
    u1[i] = r[i];
}

void mg_free(struct mg_data *d) {
  if (d != NULL) {
    struct mg_lvl **l = d->levels;
#if 0
    for (uint i = 0; i < d->nlevels; i++) {
      if (l[i]->M != NULL)
        mat_free(l[i]->M), l[i]->M = NULL;
      if (l[i]->J != NULL)
        gs_free(l[i]->J), l[i]->J = NULL;
      if (l[i]->Q != NULL)
        gs_free(l[i]->Q), l[i]->Q = NULL;
      if (l[i] != NULL)
        free(l[i]), l[i] = NULL;
    }
#endif

    if (d->levels != NULL)
      free(d->levels), d->levels = NULL;
    if (d->level_off != NULL)
      free(d->level_off), d->level_off = NULL;
    if (d->buf != NULL)
      free(d->buf), d->buf = NULL;
    free(d), d = NULL;
  }
}
