#include <genmap-impl.h>
#include <genmap-multigrid.h>

#define MIN(a, b) ((b) < (a) ? (b) : (a))

#define GS 1
#define CSR 2
#define WEIGHTED 4
#define UNWEIGHTED 8

struct csr_entry {
  ulong r, c;
  uint proc;
};

static void find_neighbors(struct array *arr, struct rsb_element *elems,
                           sint nelt, int nv, struct comm *cc, buffer *buf) {
  slong out[2][1], bfr[2][1];
  slong lelt = nelt;
  comm_scan(out, cc, gs_long, gs_add, &lelt, 1, bfr);
  ulong elem_id = out[0][0] + 1;
  ulong sequence_id = elem_id * nv;

  struct array vertices;
  size_t size = nelt * nv;
  array_init(vertex, &vertices, size);

  sint i, j;
  for (i = 0; i < lelt; i++) {
    for (j = 0; j < nv; j++) {
      vertex vrt = {.sequenceId = sequence_id,
                    .nNeighbors = 0,
                    .elementId = elem_id,
                    .vertexId = elems[i].vertices[j],
                    .workProc = elems[i].vertices[j] % cc->np};
      array_cat(vertex, &vertices, &vrt, 1);
      sequence_id++;
    }
    elem_id++;
  }

  struct crystal cr;
  crystal_init(&cr, cc);

  sarray_transfer(vertex, &vertices, workProc, 1, &cr);
  size = vertices.n;
  vertex *vPtr = vertices.ptr;

  sarray_sort(vertex, vPtr, size, vertexId, 1, buf);

  /* FIXME: Assumes quads or hexes */
  array_init(struct csr_entry, arr, 10);
  sint s = 0, e;
  struct csr_entry t;
  while (s < size) {
    e = s + 1;
    while (e < size && vPtr[s].vertexId == vPtr[e].vertexId)
      e++;
    int nnbrs = MIN(e, size) - s;

    for (i = s; i < MIN(e, size); i++) {
      t.r = vPtr[i].elementId;
      t.proc = vPtr[i].workProc;
      for (j = 0; j < nnbrs; j++) {
        t.c = vPtr[s + j].elementId;
        if (t.r != t.c)
          array_cat(struct csr_entry, arr, &t, 1);
      }
    }
    s = e;
  }

  sarray_transfer(struct csr_entry, arr, proc, 1, &cr);
  crystal_free(&cr);
  array_free(&vertices);
}

int laplacian_init(struct laplacian *l, struct rsb_element *elems, uint lelt,
                   int nv, struct comm *c, buffer *buf) {
  l->type = 0 | CSR | UNWEIGHTED;

  struct array entries;
  find_neighbors(&entries, elems, lelt, nv, c, buf);
  sarray_sort_2(struct csr_entry, entries.ptr, entries.n, r, 1, c, 1, buf);

  l->lelt = lelt;
  l->nv = nv;

  uint i, rn = 0;
  struct csr_entry *ptr = entries.ptr;
  for (i = 0; i < entries.n; i++) {}
  l->col_ids = tcalloc(ulong, entries.n);

  l->adj_off = tcalloc(uint, lelt + 1);
  l->adj_ind = tcalloc(uint, entries.n);

  l->diag = tcalloc(GenmapScalar, lelt);

  array_free(&entries);

  return 0;
}

int laplacian_weighted_init(struct laplacian *wl, struct rsb_element *elems,
                            uint lelt, int nv, struct comm *c, buffer *buf) {
  wl->type = 0 | GS | WEIGHTED;

  uint npts = nv * lelt;
  slong *vertices = tcalloc(slong, npts);
  uint i, j;
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      vertices[i * nv + j] = elems[i].vertices[j];

  wl->u = tcalloc(GenmapScalar, npts);
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      wl->u[nv * i + j] = 1.0;

  wl->gsh = gs_setup(vertices, npts, c, 0, gs_crystal_router, 0);
  gs(wl->u, gs_double, gs_add, 0, wl->gsh, buf);

  wl->diag = tcalloc(GenmapScalar, lelt);
  for (i = 0; i < lelt; i++) {
    wl->diag[i] = 0.0;
    for (j = 0; j < nv; j++)
      wl->diag[i] += wl->u[nv * i + j];
  }

  wl->lelt = lelt;
  wl->nv = nv;
  wl->col_ids = NULL;
  wl->adj_off = NULL;
  wl->adj_ind = NULL;

  if (vertices != NULL)
    free(vertices);

  return 0;
}

int laplacian_weighted(GenmapScalar *v, struct laplacian *wl, GenmapScalar *u,
                       buffer *buf) {
  uint lelt = wl->lelt;
  int nv = wl->nv;

  GenmapInt i, j;
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      wl->u[nv * i + j] = u[i];

  gs(wl->u, gs_double, gs_add, 0, wl->gsh, buf);

  for (i = 0; i < lelt; i++) {
    v[i] = wl->diag[i] * u[i];
    for (j = 0; j < nv; j++)
      v[i] -= wl->u[nv * i + j];
  }

  return 0;
}

int laplacian(GenmapScalar *v, struct laplacian *l, GenmapScalar *u,
              buffer *buf) {
  // csr_mat_apply(v, M, u, buf);
  return 0;
}

void laplacian_free(struct laplacian *l) {
  if (l->u != NULL)
    free(l->u);
  if (l->diag != NULL)
    free(l->diag);
  if (l->gsh != NULL)
    gs_free(l->gsh);
  if (l->col_ids != NULL)
    free(l->col_ids);
  if (l->adj_off != NULL)
    free(l->adj_off);
  if (l->adj_ind != NULL)
    free(l->adj_ind);
}

#undef GS
#undef CSR
#undef WEIGHTED
#undef UNWEIGHTED

#undef MIN
