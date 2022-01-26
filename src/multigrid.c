#include "multigrid.h"

#define GETPTR(p, i, off) ((char *)(p) + (off) + (i) * sizeof(entry))

struct mg_lvl {
  int nsmooth;
  GenmapScalar sigma;
  struct gs_data *J; // interpolation from level i to i+1
  struct gs_data *Q; // global to local conversion of a vector
  struct csr_laplacian *M;
};

struct mg_data {
  struct comm c;
  struct gs_data *top;
  buffer bfr;
  int nlevels;
  struct mg_lvl **levels;
  uint *level_off;
  GenmapScalar *y, *x, *b, *u, *rhs, *buf;
};

typedef struct {
  ulong r, c, rn, cn;
  uint p;
  GenmapScalar v;
} entry;

static void set_owner(char *ptr, sint n, size_t ioff, size_t ooff, slong lelg,
                      sint np) {
  sint lelt = lelg / np;
  sint nrem = lelg % np;

  ulong *iptr;
  sint *optr;
  sint i;
  slong row;
  for (i = 0; i < n; i++) {
    iptr = (ulong *)GETPTR(ptr, i, ioff);
    optr = (sint *)GETPTR(ptr, i, ooff);
    row = *iptr - 1;
    // FIXME: Assumes the 'reverse-Nek' element distribution
#if 0
    if(row<lelt*(np-nrem)) *optr=(sint) row/lelt;
    else *optr=np-nrem+(sint) (row-lelt*(np-nrem))/(lelt+1);
#else
    if (nrem == 0)
      *optr = (sint)row / lelt;
    else if (row < (lelt + 1) * nrem)
      *optr = (sint)row / (lelt + 1);
    else
      *optr = nrem + (sint)(row - (lelt + 1) * nrem) / lelt;
#endif
  }
}

static void compress_col(struct array *entries) {
  GenmapScalar v;
  sint i, j;

  i = 0;
  entry *ptr = entries->ptr;
  while (i < entries->n) {
    v = ptr[i].v, j = i + 1;
    while (j < entries->n && ptr[j].r == ptr[i].r && ptr[j].cn == ptr[i].cn)
      v += ptr[j].v, ptr[j].c = 0, j++;
    ptr[i].v = v;
    i = j;
  }
}

static void mg_lvl_setup(struct mg_data *d, uint lvl, int factor) {
  assert(lvl > 0);
  struct csr_laplacian *M0 = d->levels[lvl - 1]->M;
  uint rn0 = M0->rn, nnz0 = M0->roff[rn0];

  struct array entries = null_array;
  array_init(entry, &entries, nnz0);
  entries.n = nnz0;

  uint i, j, nn = 0;
  entry *ptr = entries.ptr;
  for (i = 0; i < rn0; i++) {
    for (j = M0->roff[i]; j < M0->roff[i + 1]; j++) {
      ptr[nn].r = M0->rstart + i;
      ptr[nn].c = M0->col[j];
      ptr[nn].rn = (ptr[nn].r + factor - 1) / factor;
      ptr[nn].cn = (ptr[nn].c + factor - 1) / factor;
      ptr[nn].v = M0->v[j];
      nn++;
    }
  }
  assert(nn == nnz0);

  struct comm *c = &d->c;
  slong out[2][1], bf[2][1];
  slong in = rn0;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bf);
  slong ng = out[1][0];

  ulong ngc = (ng / factor) > 0 ? (ng / factor) : 1;
  ulong ng_ = ngc * factor;
  if (ng > factor - 1 && ng_ < ng)
    for (j = 0; j < M0->roff[rn0]; j++) {
      if (ptr[j].c > ng_)
        ptr[j].cn -= 1;
      if (ptr[j].r > ng_)
        ptr[j].rn -= 1;
    }

  /* setup gs ids for fine level (rhs interpolation) */
  ptr = entries.ptr;
  slong *ids;
  GenmapMalloc(rn0, &ids);
  for (i = j = 0; i < nnz0; i++)
    if (ptr[i].r == ptr[i].c)
      ids[j++] = -ptr[i].cn;
  assert(j == rn0);

  /* coarsen the cols */
  buffer buf;
  buffer_init(&buf, 1024);
  if (entries.n) {
    sarray_sort_2(entry, entries.ptr, entries.n, r, 1, cn, 1, &buf);
    compress_col(&entries);
  }

  struct crystal cr;
  crystal_init(&cr, c);

  uint npc = ngc < c->np ? ngc : c->np;
  set_owner(entries.ptr, nnz0, offsetof(entry, rn), offsetof(entry, p), ngc,
            npc);
  sarray_transfer(entry, &entries, p, 1, &cr);

  // sort by rn and cn
  sarray_sort_2(entry, entries.ptr, entries.n, rn, 1, cn, 1, &buf);

  i = j = nn = 0;
  ptr = entries.ptr;
  while (i < entries.n) {
    while (j < entries.n && ptr[j].rn == ptr[i].rn)
      j++;
    i = j, nn++;
  }

  /* create the matrix */
  GenmapMalloc(1, &d->levels[lvl]);
  struct mg_lvl *l = d->levels[lvl];

  GenmapMalloc(1, &l->M);
  struct csr_laplacian *M1 = l->M;
  M1->rn = nn;
  GenmapMalloc(M1->rn + 1, &M1->roff);

  slong cn = nn;
  comm_scan(out, &d->c, gs_long, gs_add, &cn, 1, bf);
  M1->rstart = out[0][0] + 1;

  uint nnz1 = 0;
  i = j = 0;
  ptr = entries.ptr;
  while (i < entries.n) {
    while (j < entries.n && ptr[j].rn == ptr[i].rn && ptr[j].cn == ptr[i].cn)
      j++;
    i = j, nnz1++;
  }

  if (nnz1 == 0) {
    M1->col = NULL;
    M1->v = NULL;
    M1->diag = NULL;
    M1->buf = NULL;
  } else {
    GenmapMalloc(nnz1, &M1->col);
    GenmapMalloc(nnz1, &M1->v);
    GenmapMalloc(nnz1, &M1->buf);
    GenmapMalloc(M1->rn, &M1->diag);
  }

  uint rn1;
  GenmapScalar v;
  M1->roff[0] = i = j = nn = rn1 = 0;
  ptr = entries.ptr;
  while (i < entries.n) {
    v = 0.0;
    while (j < entries.n && ptr[j].rn == ptr[i].rn && ptr[j].cn == ptr[i].cn) {
      if (ptr[j].c > 0)
        v += ptr[j].v;
      j++;
    }
    M1->col[nn] = ptr[i].cn, M1->v[nn] = v, nn++;

    if ((j < entries.n && ptr[j].rn != ptr[i].rn) || j >= entries.n)
      M1->roff[++rn1] = nn;
    i = j;
  }
  assert(nn == nnz1);    // sanity check
  assert(rn1 == M1->rn); // sanity check

  /* setup gs ids for coarse level (rhs interpolation ) */
  GenmapRealloc(rn0 + rn1, &ids);
  for (i = nn = 0; i < rn1; i++)
    for (j = M1->roff[i]; j < M1->roff[i + 1]; j++)
      if (M1->rstart + i == M1->col[j]) {
        ids[rn0 + nn] = M1->col[j];
        M1->diag[i] = M1->v[j];
        nn++;
      }
  assert(nn == M1->rn);

  d->levels[lvl - 1]->J =
      gs_setup(ids, rn0 + M1->rn, &d->c, 0, gs_crystal_router, 0);

  /* setup gs handle for the mat-vec */
  GenmapRealloc(nnz1, &ids);
  for (i = 0; i < M1->rn; i++)
    for (j = M1->roff[i]; j < M1->roff[i + 1]; j++)
      if (M1->rstart + i == M1->col[j])
        ids[j] = M1->col[j];
      else
        ids[j] = -M1->col[j];

  M1->gsh = gs_setup(ids, nnz1, &d->c, 0, gs_crystal_router, 0);

  GenmapFree(ids);
  buffer_free(&buf);
  crystal_free(&cr);
  array_free(&entries);
}

struct mg_data* mg_setup(int factor, struct comm *c, struct csr_laplacian *M) {
  struct mg_data *d = tcalloc(struct mg_data, 1);
  comm_dup(&d->c, c);

  uint np = c->np, rn = M->rn;
  slong out[2][1], bf[2][1];
  slong in = rn;
  comm_scan(out, &d->c, gs_long, gs_add, &in, 1, bf);
  slong rg = out[1][0];

  d->nlevels = logbll(rg, factor) + 1;
  GenmapMalloc(d->nlevels, &d->levels);
  GenmapMalloc(d->nlevels + 1, &d->level_off);

#if 0
  if (c->id == 0)
    printf("rg = %lld factor = %d nlevels = %d\n", rg, factor, d->nlevels);
#endif

  GenmapMalloc(1, &d->levels[0]);
  d->levels[0]->M = M;
  d->levels[0]->nsmooth = 2;
  d->levels[0]->sigma = 0.6;
  d->level_off[0] = 0;
  d->level_off[1] = M->rn;

  uint i;
  uint nnz = M->roff[M->rn];
  for (i = 1; i < d->nlevels; i++) {
    mg_lvl_setup(d, i, factor);
    struct csr_laplacian *Mi = d->levels[i]->M;
    if (Mi->roff[Mi->rn] > nnz)
      nnz = Mi->roff[Mi->rn];
    d->level_off[i + 1] = d->level_off[i] + Mi->rn;
    d->levels[i]->nsmooth = 2;
    d->levels[i]->sigma = 0.6;
  }

  GenmapMalloc(d->level_off[d->nlevels], &d->x);
  GenmapMalloc(d->level_off[d->nlevels], &d->y);
  GenmapMalloc(d->level_off[d->nlevels], &d->b);
  GenmapMalloc(d->level_off[d->nlevels], &d->u);
  GenmapMalloc(d->level_off[d->nlevels], &d->rhs);
  GenmapMalloc(nnz, &d->buf);

  return d;
}

void mg_vcycle(GenmapScalar *u1, GenmapScalar *rhs, struct mg_data *d) {
  GenmapScalar *s = d->x;
  GenmapScalar *Gs = d->y;
  GenmapScalar *r = d->b;
  GenmapScalar *u = d->u;

  struct mg_lvl **lvls = d->levels;
  uint *lvl_off = d->level_off;
  struct mg_lvl *l;
  struct csr_laplacian *M;

  buffer buf;
  buffer_init(&buf, 1024);

  int nsmooth, nlevels = d->nlevels, lvl;
  GenmapScalar *diag, sigma;
  uint off, n, i, j;

  for (i = 0; i < lvl_off[nlevels]; i++)
    s[i] = Gs[i] = r[i] = u[i] = 0.0;
  for (i = 0; i < lvl_off[1]; i++)
    r[i] = rhs[i];

  for (lvl = 0; lvl < nlevels - 1; lvl++) {
    off = lvl_off[lvl];
    n = lvl_off[lvl + 1] - off;
    l = lvls[lvl];
    nsmooth = l->nsmooth;
    sigma = l->sigma;
    M = l->M;
    diag = M->diag;
    assert(n == M->rn);

    // u=sigma*D*rhs
    for (j = 0; j < n; j++)
      u[off + j] = sigma * r[off + j] / diag[j];

    // G*u
    csr_mat_apply(Gs + off, M, u + off, &buf);

    // r=rhs-Gu
    for (j = 0; j < n; j++)
      r[off + j] = r[off + j] - Gs[off + j];

    for (i = 0; i < nsmooth; i++) {
      sigma = sigma + 0.066666 / nsmooth;
      // s=sigma*D*r, u=u+s
      for (j = 0; j < n; j++) {
        s[off + j] = sigma * r[off + j] / diag[j];
        u[off + j] += s[off + j];
      }

      // G*s
      csr_mat_apply(Gs + off, M, s + off, &buf);

      // r=r-Gs
      for (j = 0; j < n; j++)
        r[off + j] = r[off + j] - Gs[off + j];
    }

    // interpolate to coarser level
    gs(r + off, gs_double, gs_add, 1, l->J, &buf);
  }

  // coarsest level
  off = lvl_off[nlevels - 1];
  n = lvl_off[nlevels] - off;

  if (n == 1) {
    l = lvls[nlevels - 1];
    M = l->M;
    assert(M->rn == 1);
    if (fabs(M->diag[0]) > sqrt(GENMAP_TOL))
      u[off] = r[off] / M->diag[0];
    else
      u[off] = 0.0;
    r[off] = u[off];
  }

  GenmapScalar over = 1.33333;
  for (lvl = nlevels - 2; lvl >= 0; lvl--) {
    l = lvls[lvl];
    off = lvl_off[lvl];
    // J*e
    gs(r + off, gs_double, gs_add, 0, l->J, &buf);

    // u = u + over * J * e
    n = lvl_off[lvl + 1] - off;
    for (j = 0; j < n; j++)
      r[off + j] = over * r[off + j] + u[off + j];
  }

  // avoid this
  for (i = 0; i < lvl_off[1]; i++)
    u1[i] = r[i];

  buffer_free(&buf);
}

void mg_free(struct mg_data *d) {
  if (d != NULL) {
    struct mg_lvl **l = d->levels;
    for (uint i = 0; i < d->nlevels; i++) {
      if (i > 0)
        csr_mat_free(l[i]->M);
      if (i < d->nlevels - 1)
        gs_free(l[i]->J);
      GenmapFree(l[i]);
    }

    // comm_free
    comm_free(&d->c);

    GenmapFree(l);
    GenmapFree(d->level_off);
    GenmapFree(d->y);
    GenmapFree(d->x);
    GenmapFree(d->b);
    GenmapFree(d->buf);
    GenmapFree(d->rhs);
    GenmapFree(d->u);

    free(d);
  }
}

#undef GETPTR
