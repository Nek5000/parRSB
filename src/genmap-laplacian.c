#include <genmap-impl.h>

#define min(a, b) ((b) < (a) ? (b) : (a))

struct array *GenmapFindNeighbors(genmap_handle h, genmap_comm c) {
  struct comm cc = c->gsc;

  sint lelt = GenmapGetNLocalElements(h);
  sint nv = GenmapGetNVertices(h);

  genmap_scan(h, c);
  ulong elem_id = GenmapGetLocalStartIndex(h) + 1;
  ulong sequenceId = elem_id * nv;

  size_t size = lelt * nv;
  struct array vertices;
  array_init(vertex, &vertices, size);

  GenmapElements elems = GenmapGetElements(h);
  sint i, j;
  for (i = 0; i < lelt; i++) {
    for (j = 0; j < nv; j++) {
      vertex vrt = {.sequenceId = sequenceId,
                    .nNeighbors = 0,
                    .elementId = elem_id,
                    .vertexId = elems[i].vertices[j],
                    .workProc = elems[i].vertices[j] % cc.np};
      array_cat(vertex, &vertices, &vrt, 1);
      sequenceId++;
    }
    elem_id++;
  }
  assert(vertices.n == lelt * nv);

  struct crystal cr;
  crystal_init(&cr, &cc);

  sarray_transfer(vertex, &vertices, workProc, 1, &cr);
  size = vertices.n;
  vertex *vPtr = vertices.ptr;

  buffer buf;
  buffer_init(&buf, 1024);
  sarray_sort(vertex, vPtr, size, vertexId, 1, &buf);

  struct array a;
  array_init(csr_entry, &a, 10);

  // FIXME: Assumes quads or hexes
  sint s = 0, e;
  csr_entry t;
  while (s < size) {
    e = s + 1;
    while (e < size && vPtr[s].vertexId == vPtr[e].vertexId)
      e++;
    int n_neighbors = min(e, size) - s;

    for (i = s; i < min(e, size); i++) {
      t.r = vPtr[i].elementId;
      t.proc = vPtr[i].workProc;
      for (j = 0; j < n_neighbors; j++) {
        t.c = vPtr[s + j].elementId;
        array_cat(csr_entry, &a, &t, 1);
      }
    }
    s = e;
  }

  sarray_transfer(csr_entry, &a, proc, 1, &cr);
  sarray_sort_2(csr_entry, a.ptr, a.n, r, 1, c, 1, &buf);
  sarray_sort(csr_entry, a.ptr, a.n, r, 1, &buf);


  struct array *nbrs = tmalloc(struct array, 1);
  array_init(entry, nbrs, lelt);

  if (a.n == 0) {
    crystal_free(&cr);
    buffer_free(&buf);
    array_free(&vertices);
    array_free(&a);
    return nbrs;
  }

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

  sarray_sort_2(entry, nbrs->ptr, nbrs->n, r, 1, c, 1, &buf);

  crystal_free(&cr);
  buffer_free(&buf);
  array_free(&vertices);
  array_free(&a);

  return nbrs;
}

int GenmapInitLaplacian(genmap_handle h, genmap_comm c) {
  struct array *entries = GenmapFindNeighbors(h, c);
  csr_mat_setup(entries, &c->gsc, &c->M);
  array_free(entries);
  free(entries);

  c->gsh = get_csr_top(c->M, &c->gsc);
  GenmapRealloc(c->M->row_off[c->M->rn], &h->b);

#if defined(GENMAP_DEBUG)
  int nnz = c->M->row_off[c->M->rn];
  double fro[2] = {0.0, 0.0}, buf[2];
  for (int i = 0; i < nnz; i++) {
    fro[0] += c->M->v[i];
    fro[1] += c->M->v[i] * c->M->v[i];
  }
  comm_allreduce(&c->gsc, gs_double, gs_add, &fro, 2, &buf);
  if (c->gsc.id == 0)
    printf("nrom(G,'1')=%g\nnorm(G,'fro')=%g\n", fro[0], fro[1]);
#endif

  return 0;
}

int GenmapLaplacian(genmap_handle h, genmap_comm c, GenmapScalar *u,
                    GenmapScalar *v) {
  csr_mat_gather(c->M, c->gsh, u, h->b, &h->buf);
  csr_mat_apply(v, c->M, h->b);

  return 0;
}

#undef min
