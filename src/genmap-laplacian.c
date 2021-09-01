#include <genmap-multigrid-precon.h>
#include <genmap-impl.h>

#define MIN(a, b) ((b) < (a) ? (b) : (a))

static void genmap_find_neighbors(struct array *nbrs, struct rsb_element *elems,
                                  sint nelt, int nv, struct comm *cc,
                                  buffer *buf) {
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

int GenmapInitLaplacian(struct csr_mat *M, struct rsb_element *elems, uint lelt,
                        int nv, struct comm *c, buffer *buf) {
  struct array entries;
  genmap_find_neighbors(&entries, elems, lelt, nv, c, buf);
  csr_mat_setup(M, &entries, c, buf);
  array_free(&entries);

  return 0;
}

int GenmapLaplacian(GenmapScalar *v, struct csr_mat *M, GenmapScalar *u,
                    buffer *buf) {
  csr_mat_apply(v, M, u, buf);
  return 0;
}

int GenmapInitLaplacianWeighted(struct laplacian *gl, struct rsb_element *elems,
                                uint lelt, int nv, struct comm *c, buffer *buf) {
  uint npts = nv * lelt;

  slong *vertices = tcalloc(slong, npts);
  GenmapInt i, j;
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      vertices[i * nv + j] = elems[i].vertices[j];

  gl->u = tcalloc(GenmapScalar, npts);
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      gl->u[nv * i + j] = 1.0;

  gl->gsh = gs_setup(vertices, npts, c, 0, gs_crystal_router, 0);
  gs(gl->u, gs_double, gs_add, 0, gl->gsh, buf);

  gl->weights = tcalloc(GenmapScalar, lelt);
  for (i = 0; i < lelt; i++) {
    gl->weights[i] = 0.0;
    for (j = 0; j < nv; j++)
      gl->weights[i] += gl->u[nv * i + j];
  }

  gl->lelt = lelt;
  gl->nv = nv;

  if (vertices != NULL)
    free(vertices);

  return 0;
}

int GenmapLaplacianWeighted(GenmapScalar *v, struct laplacian *gl,
                            GenmapScalar *u, buffer *buf) {
  uint lelt = gl->lelt;
  int nv = gl->nv;

  GenmapInt i, j;
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      gl->u[nv * i + j] = u[i];

  gs(gl->u, gs_double, gs_add, 0, gl->gsh, buf);

  for (i = 0; i < lelt; i++) {
    v[i] = gl->weights[i] * u[i];
    for (j = 0; j < nv; j++)
      v[i] -= gl->u[nv * i + j];
  }

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

#undef MIN
