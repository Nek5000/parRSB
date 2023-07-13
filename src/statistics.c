#include "parrsb-impl.h"

static uint get_partition(const struct comm *gc, const struct comm *lc) {
  // Find the partition id. A partition is a group of processors sharing the
  // same local communicator.
  sint out[2][1], wrk[2][1], root = (lc->id == 0);
  comm_scan(out, gc, gs_int, gs_add, &root, 1, wrk);

  sint smin = out[0][0], smax = out[0][0];
  comm_allreduce(lc, gs_int, gs_min, &smin, 1, wrk);
  comm_allreduce(lc, gs_int, gs_max, &smax, 1, wrk);
  assert(smin == smax && "MPI processes in a partition must have the same id.");

  return smin;
}

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

  const struct rsb_element *const pe =
      (const struct rsb_element *const)elems->ptr;
  struct vertex_t vt = {.partition = get_partition(gc, lc)};
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

static struct array pgeom;
static struct comm comm;
static uint pgeom_initialized = 0;
static uint ndim = 0;

struct pgeom_t {
  uint partition, level;
  double centroid[3], min[3], max[3];
  uint p;
};

void dump_part_start(const struct comm *gc, const uint nv) {
  assert(!pgeom_initialized && "Partition geometry is already initialized.");

  comm_dup(&comm, gc);
  ndim = (nv == 8) ? 3 : 2;
  array_init(struct pgeom_t, &pgeom, 1024);
  pgeom_initialized = 1;
}

void dump_part_geom(const struct comm *lc, const struct array *const elems,
                    const int level, buffer *bfr) {
  assert(pgeom_initialized && "Partition geometry is not initialized.");

  const struct rsb_element *const pe =
      (const struct rsb_element *const)elems->ptr;

  // Find the centroid and the bounding box of the partition.
  double centroid[3] = {0.0, 0.0, 0.0};
  double max[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
  double min[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
  const uint n = elems->n;
  for (uint e = 0; e < n; e++) {
    for (uint d = 0; d < ndim; d++) {
      double c = pe[e].coord[d];
      centroid[d] += c;
      max[d] = (max[d] < c) ? c : max[d];
      min[d] = (min[d] > c) ? c : min[d];
    }
  }
  for (uint d = 0; d < ndim; d++)
    centroid[d] /= n;

  double wrk[3];
  comm_allreduce(lc, gs_double, gs_min, min, ndim, wrk);
  comm_allreduce(lc, gs_double, gs_max, max, ndim, wrk);
  comm_allreduce(lc, gs_double, gs_add, centroid, ndim, wrk);
  for (uint d = 0; d < ndim; d++)
    centroid[d] /= lc->np;

  // Partition root accumulates the partition geometry.
  if (lc->id == 0) {
    struct pgeom_t pg = {.partition = get_partition(&comm, lc),
                         .level = level,
                         .centroid = {centroid[0], centroid[1], centroid[2]},
                         .max = {max[0], max[1], max[2]},
                         .min = {min[0], min[1], min[2]},
                         .p = 0};
    array_cat(struct pgeom_t, &pgeom, &pg, 1);
  }
}

void dump_part_end(const char *prefix) {
  assert(pgeom_initialized && "Partition geometry is not initialized.");

  const uint size = strnlen(prefix, 64);
  assert(size < 64 && "Prefix must be less than 64 characters.");

  // Send all the data to global root.
  struct crystal cr;
  crystal_init(&cr, &comm);
  sarray_transfer(struct pgeom_t, &pgeom, p, 0, &cr);
  crystal_free(&cr);

  if (comm.id == 0) {
    const char name[BUFSIZ];
    snprintf((char *)name, BUFSIZ, "%s_partition_geom.txt", prefix);

    FILE *fp = fopen(name, "w");
    if (!fp) {
      fprintf(stderr, "Failed to open %s for writing.\n", name);
      exit(EXIT_FAILURE);
    }

    buffer bfr;
    buffer_init(&bfr, 1024);
    sarray_sort_2(struct pgeom_t, pgeom.ptr, pgeom.n, level, 0, partition, 0,
                  &bfr);
    buffer_free(&bfr);

    const struct pgeom_t *const pg = (const struct pgeom_t *const)pgeom.ptr;
    for (uint i = 0; i < pgeom.n; i++) {
      fprintf(fp, "%u %u %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", pg[i].level,
              pg[i].partition, pg[i].centroid[0], pg[i].centroid[1],
              pg[i].centroid[2], pg[i].min[0], pg[i].min[1], pg[i].min[2],
              pg[i].max[0], pg[i].max[1], pg[i].max[2]);
    }
    fclose(fp);
  }

  comm_free(&comm);
  array_free(&pgeom);
  ndim = pgeom_initialized = 0;
}
