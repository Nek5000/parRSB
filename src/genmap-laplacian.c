#include <genmap-impl.h>
#include <genmap-multigrid.h>

#define MIN(a, b) ((b) < (a) ? (b) : (a))

#define GS 1
#define CSR 2
#define WEIGHTED 4
#define UNWEIGHTED 8

typedef struct {
  ulong r, c;
  uint proc;
} csr_entry;

typedef struct {
  ulong r, c;
  uint p;
  GenmapScalar v;
} entry;

static void find_neighbors(struct array *nbrs, struct rsb_element *elems,
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

  struct array a;
  array_init(csr_entry, &a, 10);

  /* FIXME: Assumes quads or hexes */
  sint s = 0, e;
  csr_entry t;
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
        array_cat(csr_entry, &a, &t, 1);
      }
    }
    s = e;
  }

  sarray_transfer(csr_entry, &a, proc, 1, &cr);
  array_init(entry, nbrs, lelt);
  if (a.n == 0) {
    crystal_free(&cr);
    array_free(&vertices);
    array_free(&a);
  }

  sarray_sort_2(csr_entry, a.ptr, a.n, r, 1, c, 1, buf);
  csr_entry *aptr = a.ptr;
  entry *nptr = nbrs->ptr;

  entry ee = {0, 0, 0, 0, 0, 0.0}, ep = {0, 0, 0, 0, 0.0};
  ep.r = aptr[0].r;
  ep.c = aptr[0].c;
  array_cat(entry, nbrs, &ep, 1);

  for (i = 1; i < a.n; i++) {
    ee.r = aptr[i].r;
    ee.c = aptr[i].c;
    if (ee.r != ep.r || ee.c != ep.c) {
      array_cat(entry, nbrs, &ee, 1);
      ep = ee;
    }
  }

  sarray_sort_2(entry, nbrs->ptr, nbrs->n, r, 1, c, 1, buf);

  crystal_free(&cr);
  array_free(&vertices);
  array_free(&a);
}

int laplacian_init(struct laplacian *l, struct rsb_element *elems, uint lelt,
                   int nv, struct comm *c, buffer *buf) {
  l->type = 0 | CSR | UNWEIGHTED;

  struct array entries;
  find_neighbors(&entries, elems, lelt, nv, c, buf);

  struct csr_mat M;
  csr_mat_setup(&M, &entries, c, buf);

  l->lelt = lelt;
  l->nv = nv;

  l->off = tcalloc(uint, M.rn + 1);
  if (l->off != NULL)
    memcpy(l->off, M.roff, sizeof(uint) * (M.rn + 1));

  l->weights = tcalloc(GenmapScalar, M.rn);
  if (l->weights != NULL)
    memcpy(l->weights, M.diag, sizeof(GenmapScalar) * M.rn);

  l->u = tcalloc(GenmapScalar, M.roff[M.rn]);
  if (l->u != NULL)
    memcpy(l->u, M.v, sizeof(GenmapScalar) * M.roff[M.rn]);

  csr_mat_free(&M);
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

  wl->weights = tcalloc(GenmapScalar, lelt);
  for (i = 0; i < lelt; i++) {
    wl->weights[i] = 0.0;
    for (j = 0; j < nv; j++)
      wl->weights[i] += wl->u[nv * i + j];
  }

  wl->lelt = lelt;
  wl->nv = nv;

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
    v[i] = wl->weights[i] * u[i];
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
  if (l->weights != NULL)
    free(l->weights);
  if (l->gsh != NULL)
    gs_free(l->gsh);
}

#undef GS
#undef CSR
#undef WEIGHTED
#undef UNWEIGHTED

#undef MIN
