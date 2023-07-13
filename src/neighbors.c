#include "parrsb-impl.h"

uint get_neighbors(const struct array *const elems, const unsigned nv,
                   const struct comm *gc, const struct comm *lc, buffer *bfr) {
  const uint n = elems->n;
  const uint size = elems->n * nv;

  struct vertex_t {
    ulong v;
    uint p, partition;
  };

  struct array vertices;
  array_init(struct vertex_t, &vertices, size);

  // Find the partition id. A partition is a group of processors sharing the
  // same local communicator.
  sint out[2][1], wrk[2][1], root = (lc->id == 0);
  comm_scan(out, gc, gs_int, gs_add, &root, 1, wrk);

  sint smin = out[0][0], smax = out[0][0];
  comm_allreduce(lc, gs_int, gs_min, &smin, 1, wrk);
  comm_allreduce(lc, gs_int, gs_max, &smax, 1, wrk);
  assert(smin == smax && "MPI processes in a partition must have the same id.");

  const struct rsb_element *const pe =
      (const struct rsb_element *const)elems->ptr;
  struct vertex_t vt = {.partition = smin};
  for (uint i = 0; i < n; i++) {
    for (uint v = 0; v < nv; v++) {
      vt.v = pe[i].vertices[v], vt.p = vt.v % gc->np;
      array_cat(struct vertex_t, &vertices, &vt, 1);
    }
  }

  struct crystal cr;
  crystal_init(&cr, gc);

  sarray_transfer(struct vertex_t, &vertices, p, 1, &cr);
  sarray_sort(struct vertex_t, vertices.ptr, vertices.n, v, 1, bfr);

  struct array neighbors;
  array_init(struct vertex_t, &neighbors, vertices.n * 27);

  const struct vertex_t *const pv = (const struct vertex_t *const)vertices.ptr;
  uint s = 0;
  while (s < vertices.n) {
    uint e = s + 1;
    while (e < vertices.n && pv[s].v == pv[e].v)
      e++;
    for (uint i = s; i < e; i++) {
      struct vertex_t vt = pv[s];
      for (uint j = s; j < e; j++) {
        vt.partition = pv[j].partition;
        array_cat(struct vertex_t, &neighbors, &vt, 1);
      }
    }
    s = e;
  }
  array_free(&vertices);

  sarray_transfer(struct vertex_t, &neighbors, p, 0, &cr);
  crystal_free(&cr);
  sarray_sort(struct vertex_t, neighbors.ptr, neighbors.n, partition, 0, bfr);

  uint nneighbors = 0;
  if (neighbors.n > 0) {
    nneighbors = 1;
    const struct vertex_t *const pn =
        (const struct vertex_t *const)neighbors.ptr;
    for (uint i = 1; i < neighbors.n; i++) {
      if (pn[i].partition > pn[i - 1].partition)
        nneighbors++;
    }
  }

  array_free(&neighbors);
  return nneighbors;
}
