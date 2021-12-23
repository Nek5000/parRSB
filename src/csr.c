#include "genmap-impl.h"
#include "multigrid.h"

#define ABS(i) ((i < 0) ? -i : i)

typedef struct {
  ulong r, c;
  uint p;
  GenmapScalar v;
} entry;

static struct gs_data *get_csr_top(struct csr_laplacian *M, struct comm *c) {
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

void csr_mat_gather(GenmapScalar *buf, struct csr_laplacian *M, GenmapScalar *x,
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

void csr_mat_apply(GenmapScalar *y, struct csr_laplacian *M, GenmapScalar *x,
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

int csr_mat_free(struct csr_laplacian *M) {
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
