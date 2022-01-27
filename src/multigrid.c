#include "multigrid.h"
#include <math.h>

struct mg_lvl {
  int nsmooth;
  scalar sigma;
  struct gs_data *J; // Interpolation from level i to i + 1
  struct gs_data *Q; // gs handle for matrix vector product
  struct mat *M;
};

struct mg_data {
  int nlevels;
  struct mg_lvl **levels;
  uint *level_off;
  scalar *buf;
};

// These definitions should match what is defined in coarse.c
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

static int logbll(slong n, int a) {
  assert(a > 0);

  int k = 0;
  while (n > 1)
    n /= a, k++;

  return k;
}

static struct gs_data *setup_Q(const struct mat *M, const struct comm *c,
                               buffer *bfr) {
  // Setup gs handle for the mat-vec
  uint nnz = M->n > 0 ? M->Lp[M->n] + (M->D != NULL) * M->n : 0;
  buffer_reserve(bfr, nnz * sizeof(slong));

  slong *ids = (slong *)bfr->ptr;
  uint i, j;
  for (i = 0; i < M->n; i++)
    for (j = M->Lp[i]; j < M->Lp[i + 1]; j++)
      ids[j] = -M->col[M->Li[j]];
  for (i = 0; i < M->n; i++)
    ids[j++] = M->start + i;

  return gs_setup(ids, nnz, c, 0, gs_crystal_router, 0);
}

static void mg_lvl_setup(struct mg_data *d, const uint lvl, const int factor,
                         const struct comm *c, struct crystal *cr,
                         struct array *entries, buffer *bfr) {
  assert(lvl > 0);
  struct mat *M = d->levels[lvl - 1]->M;
  uint nnz = M->n > 0 ? M->Lp[M->n] + (M->D != NULL) * M->n : 0;

  // Reserve enough space in the array
  array_reserve(struct mat_ij, entries, nnz);
  entries->n = 0;

  // Reserve enough memory for gs ids for interpolation
  slong *ids = (slong *)tcalloc(slong, 2 * M->n);

  // Calculate coarse level parameters: ngc, npc, nelt, nrem
  slong out[2][1], bf[2][1], in = M->n;
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
  for (i = 0; i < M->n; i++) {
    ulong r = M->start + i;
    for (j = M->Lp[i]; j < M->Lp[i + 1]; j++) {
      m.r = (r + factor - 1) / factor;
      if (m.r > ngc)
        m.r--;
      m.c = (M->col[M->Li[j]] + factor - 1) / factor;
      if (m.c > ngc)
        m.c--;
      m.v = M->L[j];
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
    m.v = M->D[i];
    SETP(m.p, m.r, npc, nelt, nrem);
    ids[k++] = -m.c;
    array_cat(struct mat_ij, entries, &m, 1);
  }

#undef SETP

  sarray_transfer(struct mat_ij, entries, p, 0, cr);
  sarray_sort_2(struct mat_ij, entries->ptr, entries->n, r, 1, c, 1, bfr);

  struct mg_lvl *l = d->levels[lvl] = tcalloc(struct mg_lvl, 1);
  M = l->M = csr_setup_ext(entries, 1, bfr);

  // Setup gs ids for coarse level (rhs interpolation )
  for (i = 0; i < M->n; i++)
    ids[k++] = M->start + i;

  d->levels[lvl - 1]->J = gs_setup(ids, k, c, 0, gs_crystal_router, 0);
  free(ids);
  d->levels[lvl - 1]->Q = setup_Q(M, c, bfr);
}

struct mg_data *mg_setup(const uint nelt, const ulong *eid, const slong *vtx,
                         int nv, int factor, struct comm *c, buffer *bfr) {
  slong out[2][1], bf[2][1], in = nelt;
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
#if 0
  for (uint i = 0; i < nelt; i++)
    printf("e1 = %lld\n", eid[i]);
#endif
  d->levels[0]->M = csr_from_conn(nelt, eid, vtx, nv, 1, c, bfr);

#if 0
  mat_print(d->levels[0]->M);
  fflush(stdout);
#endif
  d->levels[0]->Q = setup_Q(d->levels[0]->M, c, bfr);

  d->level_off = tcalloc(uint, d->nlevels + 1);
  d->level_off[0] = 0;
  d->level_off[1] = nelt;

  struct crystal cr;
  crystal_init(&cr, c);

  struct mat *M = d->levels[0]->M;
  uint nnz = M->n > 0 ? M->Lp[M->n] : 0;

  struct array entries;
  array_init(struct mat_ij, &entries, nnz);

  for (uint i = 1; i < d->nlevels; i++) {
    mg_lvl_setup(d, i, factor, c, &cr, &entries, bfr);
    M = d->levels[i]->M;
#if 0
    if (c->id == 0)
      printf("i = %d\n", i);
    mat_print(d->levels[i]->M);
    fflush(stdout);
#endif
    if (M->n > 0 && M->Lp[M->n] + M->n > nnz)
      nnz = M->Lp[M->n] + M->n;
    d->level_off[i + 1] = d->level_off[i] + M->n;
    d->levels[i]->nsmooth = 2;
    d->levels[i]->sigma = 0.6;
  }

  d->buf = tcalloc(scalar, 5 * d->level_off[d->nlevels] + nnz);

  array_free(&entries);
  crystal_free(&cr);

  return d;
}

static void mat_vec(scalar *y, struct mat *M, scalar *x, struct gs_data *gsh,
                    scalar *buf, buffer *bfr) {
  uint n = M->n > 0 ? M->n : 0;
  uint nnz = n > 0 ? M->Lp[n] : 0;
  uint i;
  for (i = 0; i < nnz; i++)
    buf[i] = 0.0; // Is this really necessary?
  for (i = 0; i < M->n; i++)
    y[i] = buf[nnz + i] = x[i];
  gs(buf, gs_double, gs_add, 0, gsh, bfr);

  uint j, je;
  for (i = 0; i < n; i++) {
    for (y[i] *= M->D[i], j = M->Lp[i], je = M->Lp[i + 1]; j != je; j++)
      y[i] += M->L[j] * buf[j];
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
    struct mat *M = l->M;

    // u = sigma*inv(D)*rhs
    scalar sigma = l->sigma;
    for (j = 0; j < n; j++)
      u[off + j] = sigma * r[off + j] / M->D[j];

    // G*u
    mat_vec(Gs + off, M, u + off, l->Q, buf, bfr);

    // r = rhs - Gu
    for (j = 0; j < n; j++)
      r[off + j] = r[off + j] - Gs[off + j];

    for (i = 0; i < l->nsmooth; i++) {
      sigma = sigma + 0.066666 / l->nsmooth;

      // s = sigma*inv(D)*r, u = u + s
      for (j = 0; j < n; j++) {
        s[off + j] = sigma * r[off + j] / M->D[j];
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
    struct mat *M = l->M;
    assert(M->n == 1);
    if (fabs(M->D[0]) > sqrt(TOL))
      u[off] = r[off] / M->D[0];
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
    for (uint i = 0; i < d->nlevels; i++) {
      mat_free(l[i]->M);
      if (l[i]->J != NULL)
        gs_free(l[i]->J);
      if (l[i]->Q != NULL)
        gs_free(l[i]->Q);
      free(l[i]);
    }

    free(d->levels);
    free(d->level_off);
    free(d->buf);
    free(d);
  }
}
