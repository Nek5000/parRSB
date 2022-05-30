#include "multigrid.h"
#include <math.h>

struct mg_lvl {
  int nsmooth;
  scalar sigma;
  struct gs_data *J; // Interpolation from level i to i + 1
  struct gs_data *Q; // gs handle for matrix vector product
  struct par_mat *M;
};

struct mg {
  int nlevels;
  struct mg_lvl **levels;
  uint *level_off;
  scalar *buf;
};

static int logbll(slong n, int a) {
  assert(a > 0);

  int k = (n > 0);
  while (n > 1)
    n = (n + a - 1) / a, k++;

  return k;
}

//=============================================================================
// MG setup
//
static scalar sigma_cheb(int k, int n, scalar lmin, scalar lmax) {
  k = (k - 1) % n + 1;
  scalar theta = M_PI * (k - 0.5) / n;
  scalar lamk = lmin + 0.5 * (lmax - lmin) * (cos(theta) + 1);
  return 1 / lamk;
}

static void mg_setup_aux(struct mg *d, const uint lvl, const int factor,
                         const struct comm *c, struct crystal *cr,
                         struct array *entries, buffer *bfr) {
  assert(lvl > 0);
  struct par_mat *M = d->levels[lvl - 1]->M;
  uint nnz = M->rn > 0 ? M->adj_off[M->rn] + IS_DIAG(M) * M->rn : 0;

  // Reserve enough space in the array
  array_reserve(struct mij, entries, nnz);
  entries->n = 0;

  // Reserve enough memory for gs ids for interpolation
  slong *ids = (slong *)tcalloc(slong, 2 * M->rn);

  // Calculate the minimum row id
  slong rs = M->rn > 0 ? M->rows[0] : LLONG_MAX, b;
  comm_allreduce(c, gs_long, gs_min, &rs, 1, &b);
  rs -= 1;

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
        p = s + (r - t - 1) / (nelt + 1);                                      \
    }                                                                          \
    assert(p >= 0 && p < np);                                                  \
  }

  uint i, j, k = 0, n = 0;
  struct mij m;
  for (i = 0; i < M->rn; i++) {
    ulong r = M->rows[i] - rs;
    for (j = M->adj_off[i]; j < M->adj_off[i + 1]; j++) {
      m.r = (r + factor - 1) / factor;
      if (m.r > ngc)
        m.r--;
      m.c = (M->cols[M->adj_idx[j]] - rs + factor - 1) / factor;
      if (m.c > ngc)
        m.c--;
      m.v = M->adj_val[j];
      assert(m.r > 0);
      assert(m.c > 0);
      SETP(m.p, m.r, npc, nelt, nrem);
      array_cat(struct mij, entries, &m, 1);
    }
    m.r = (r + factor - 1) / factor;
    if (m.r > ngc)
      m.r--;
    m.c = m.r;
    assert(m.r > 0);
    m.v = M->diag_val[i];
    SETP(m.p, m.r, npc, nelt, nrem);
    ids[k++] = -m.c;
    array_cat(struct mij, entries, &m, 1);
  }

#undef SETP

  sarray_transfer(struct mij, entries, p, 0, cr);
  sarray_sort_2(struct mij, entries->ptr, entries->n, r, 1, c, 1, bfr);

  struct mg_lvl *l = d->levels[lvl] = tcalloc(struct mg_lvl, 1);
  M = l->M = par_csr_setup_ext(entries, 1, bfr);

  // Setup gs ids for coarse level (rhs interpolation )
  ids = (slong *)trealloc(slong, ids, k + M->rn);
  for (i = 0; i < M->rn; i++)
    ids[k++] = M->rows[i];
  d->levels[lvl - 1]->J = gs_setup(ids, k, c, 0, gs_pairwise, 0);
  free(ids);

  d->levels[lvl]->Q = setup_Q(M, c, bfr);
}

struct mg *mg_setup(const struct par_mat *M, const int factor,
                    const struct comm *c, struct crystal *cr, buffer *bfr) {
  assert(IS_CSR(M));
  assert(M->rn == 0 || IS_DIAG(M));

  slong out[2][1], wrk[2][1], in = M->rn;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong rg = out[1][0];

  struct mg *d = (struct mg *)tcalloc(struct mg, 1);
  d->nlevels = logbll(rg, factor);

  if (d->nlevels == 0) {
    d->levels = NULL;
    d->level_off = NULL;
    d->buf = NULL;
    return d;
  }

  d->levels = (struct mg_lvl **)tcalloc(struct mg_lvl *, d->nlevels);
  d->level_off = tcalloc(uint, d->nlevels + 1);

  // Setup Level 1, keeps a pointer to input matrix
  d->levels[0] = (struct mg_lvl *)tcalloc(struct mg_lvl, 1);
  d->levels[0]->nsmooth = 3;
  d->levels[0]->sigma = sigma_cheb(1, d->levels[0]->nsmooth + 1, 1, 2);
  d->levels[0]->M = (struct par_mat *)M;
  d->levels[0]->Q = setup_Q(M, c, bfr);

  d->level_off[0] = 0;
  d->level_off[1] = M->rn;

  uint nnz = M->rn > 0 ? M->adj_off[M->rn] + M->rn : 0;
  struct array mijs;
  array_init(struct mij, &mijs, nnz);

  uint l = 1;
  for (; l < d->nlevels; l++) {
    mg_setup_aux(d, l, factor, c, cr, &mijs, bfr);
    struct par_mat *Ml = d->levels[l]->M;
    if (Ml->rn > 0 && Ml->adj_off[Ml->rn] + Ml->rn > nnz)
      nnz = Ml->adj_off[Ml->rn] + Ml->rn;
    d->level_off[l + 1] = d->level_off[l] + Ml->rn;
    d->levels[l]->nsmooth = 3;
    d->levels[l]->sigma = sigma_cheb(1, d->levels[l]->nsmooth + 1, 1, 2);
  }

  d->levels[l - 1]->J = NULL;
  d->buf = tcalloc(scalar, 5 * d->level_off[d->nlevels] + nnz);

  array_free(&mijs);

  return d;
}

//==============================================================================
// MG V-cycle and related functions
//
void mg_vcycle(scalar *u1, scalar *rhs, struct mg *d, struct comm *c,
               buffer *bfr) {
  if (d == NULL)
    return;

  int nlevels = d->nlevels;
  uint *lvl_off = d->level_off, nnz = lvl_off[nlevels], i;
  scalar *s = d->buf, *Gs = s + nnz, *r = Gs + nnz, *u = r + nnz;
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
    assert(IS_CSR(M));
    assert(n == 0 || IS_DIAG(M));
    assert(n == M->rn);

    // u = sigma*inv(D)*rhs
    scalar sigma = l->sigma;
    for (j = 0; j < n; j++)
      u[off + j] = sigma * r[off + j] / M->diag_val[j];

    // G*u
    mat_vec_csr(Gs + off, u + off, M, l->Q, buf, bfr);

    // r = rhs - Gu
    for (j = 0; j < n; j++)
      r[off + j] = r[off + j] - Gs[off + j];

    for (i = 0; i < l->nsmooth - 1; i++) {
      sigma = sigma + 0.066666 / l->nsmooth;

      // s = sigma*inv(D)*r, u = u + s
      for (j = 0; j < n; j++) {
        s[off + j] = sigma * r[off + j] / M->diag_val[j];
        u[off + j] += s[off + j];
      }

      // G*s
      mat_vec_csr(Gs + off, s + off, M, l->Q, buf, bfr);

      // r = r - Gs
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
    if (fabs(M->diag_val[0]) > 1e-6)
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

void mg_free(struct mg *d) {
  if (d != NULL) {
    struct mg_lvl **l = d->levels;
    for (uint i = 0; i < d->nlevels; i++) {
      if (l[i]->M != NULL)
        par_mat_free(l[i]->M), l[i]->M = NULL;
      if (l[i]->J != NULL)
        gs_free(l[i]->J), l[i]->J = NULL;
      if (l[i]->Q != NULL)
        gs_free(l[i]->Q), l[i]->Q = NULL;
      if (l[i] != NULL)
        free(l[i]), l[i] = NULL;
    }

    if (d->levels != NULL)
      free(d->levels), d->levels = NULL;
    if (d->level_off != NULL)
      free(d->level_off), d->level_off = NULL;
    if (d->buf != NULL)
      free(d->buf), d->buf = NULL;
    // We don't set d to NULL here -- we need to pass struct `struct mg **d` to
    // mg_free in order to do so
    free(d);
  }
}
