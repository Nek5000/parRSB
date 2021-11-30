#include <genmap-impl.h>
#include <genmap-multigrid.h>

#define ABS(i) ((i < 0) ? -i : i)

typedef struct {
  ulong r, c;
  uint p;
  GenmapScalar v;
} entry;

static struct gs_data *get_csr_top(struct csr_mat *M, struct comm *c) {
  const uint rn = M->rn;
  const uint n = M->roff[rn];

  slong *ids;
  if (n > 0)
    GenmapMalloc(n, &ids);

  uint i, j;
  for (i = 0; i < rn; i++)
    for (j = M->roff[i]; j < M->roff[i + 1]; j++)
      if (M->rstart + i == M->col[j])
        ids[j] = M->col[j];
      else
        ids[j] = -M->col[j];

  struct gs_data *gsh = gs_setup(ids, n, c, 0, gs_pairwise, 0);

  if (n > 0)
    GenmapFree(ids);

  return gsh;
}

void csr_mat_setup(struct csr_mat *M, struct array *entries, struct comm *c,
                   buffer *buf) {
  entry *ptr = entries->ptr;
  sarray_sort_2(entry, ptr, entries->n, r, 1, c, 1, buf);

  uint i = 0, j, n = 0;
  while (i < entries->n) {
    j = i + 1;
    while (j < entries->n && ptr[i].r == ptr[j].r)
      j++;
    i = j, n++;
  }

  M->rn = n;
  GenmapMalloc(n + 1, &M->roff);

  if (n == 0) {
    M->col = NULL;
    M->v = NULL;
    M->buf = NULL;
    M->diag = NULL;
  } else {
    GenmapMalloc(entries->n, &M->col);
    GenmapMalloc(entries->n, &M->v);
    GenmapMalloc(entries->n, &M->buf);
    GenmapMalloc(M->rn, &M->diag);
  }

  slong out[2][1], bf[2][1];
  slong in = M->rn;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bf);
  M->rstart = out[0][0] + 1;

  ptr = entries->ptr;
  uint rn = 0;
  for (i = 0; i < entries->n; i++) {
    if (ptr[i].r == ptr[i].c)
      M->diag[rn++] = ptr[i].v;
    M->col[i] = ptr[i].c;
    M->v[i] = ptr[i].v;
  }
  printf("rn = %d, M->rn = %d\n", rn, M->rn);
  assert(rn == M->rn);

  M->roff[0] = 0;
  i = 0;
  uint nn = 0;
  while (i < entries->n) {
    j = i + 1;
    while (j < entries->n && ptr[i].r == ptr[j].r)
      j++;
    i = M->roff[++nn] = j;
  }
  assert(n == nn);
  assert(M->roff[n] == entries->n);

  M->gsh = get_csr_top(M, c);
}

void csr_mat_gather(GenmapScalar *buf, struct csr_mat *M, GenmapScalar *x,
                    buffer *bfr) {
  ulong s = M->rstart;
  sint i, j;
  for (i = 0; i < M->rn; i++)
    for (j = M->roff[i]; j < M->roff[i + 1]; j++)
      if (M->col[j] == s + i)
        buf[j] = x[i];
      else
        buf[j] = 0.0;

  gs(buf, gs_scalar, gs_add, 0, M->gsh, bfr);
}

void csr_mat_apply(GenmapScalar *y, struct csr_mat *M, GenmapScalar *x,
                   buffer *bfr) {
  csr_mat_gather(M->buf, M, x, bfr);

  const uint rn = M->rn;
  const uint *offsets = M->roff;
  const GenmapScalar *v = M->v;

  uint i, j, je;
  for (i = 0; i < rn; i++) {
    for (y[i] = 0.0, j = offsets[i], je = offsets[i + 1]; j < je; j++)
      y[i] += M->v[j] * M->buf[j];
  }
}

void csr_mat_print(struct csr_mat *M, struct comm *c) {
  const sint rn = M->rn;
  const uint *offsets = M->roff;
  const GenmapScalar *v = M->v;
  const ulong *col = M->col;

  uint i, j, k;

  for (k = 0; k < c->np; k++) {
    genmap_barrier(c);
    if (c->id == k) {
      for (i = 0; i < rn; i++)
        fprintf(stderr, "csr %lld: %.10lf\n", M->rstart + i, M->diag[i]);
    }
    fflush(stderr);
  }
}

int csr_mat_free(struct csr_mat *M) {
  if (M->col)
    GenmapFree(M->col);
  if (M->v)
    GenmapFree(M->v);
  if (M->buf)
    GenmapFree(M->buf);
  if (M->diag)
    GenmapFree(M->diag);
  if (M->roff)
    GenmapFree(M->roff);
  if (M->gsh)
    gs_free(M->gsh);

  return 0;
}
