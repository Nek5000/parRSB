#include "metrics.h"
#include "parrsb-impl.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

parrsb_options parrsb_default_options = {
    // General options
    .partitioner = 0,
    .levels = 2,
    .repair = 0,
    .verbose_level = 0,
    .profile_level = 0,
    // RSB common (Lanczos and MG) options
    .rsb_algo = 0,
    .rsb_pre = 1,
    .rsb_max_iter = 50,
    .rsb_max_passes = 50,
    .rsb_tol = 1e-5,
    .rsb_dump_stats = 0,
    // RSB MG specific options
    .rsb_mg_grammian = 0,
    .rsb_mg_factor = 2,
    .rsb_mg_sagg = 0};

static char *ALGO[3] = {"RSB", "RCB", "RIB"};

static void update_options(parrsb_options *const options) {
#define UPDATE_OPTION(OPT, STR, IS_INT)                                        \
  do {                                                                         \
    const char *val = getenv(STR);                                             \
    if (val != NULL) {                                                         \
      if (IS_INT)                                                              \
        options->OPT = atoi(val);                                              \
      else                                                                     \
        options->OPT = atof(val);                                              \
    }                                                                          \
  } while (0)

  UPDATE_OPTION(partitioner, "PARRSB_PARTITIONER", 1);
  UPDATE_OPTION(levels, "PARRSB_LEVELS", 1);
  UPDATE_OPTION(repair, "PARRSB_REPAIR", 1);
  UPDATE_OPTION(verbose_level, "PARRSB_VERBOSE_LEVEL", 1);
  UPDATE_OPTION(profile_level, "PARRSB_PROFILE_LEVEL", 1);
  UPDATE_OPTION(rsb_algo, "PARRSB_RSB_ALGO", 1);
  UPDATE_OPTION(rsb_pre, "PARRSB_RSB_PRE", 1);
  UPDATE_OPTION(rsb_max_iter, "PARRSB_RSB_MAX_ITER", 1);
  UPDATE_OPTION(rsb_max_passes, "PARRSB_RSB_MAX_PASSES", 1);
  UPDATE_OPTION(rsb_tol, "PARRSB_RSB_TOL", 0);
  UPDATE_OPTION(rsb_dump_stats, "PARRSB_DUMP_STATS", 1);
  UPDATE_OPTION(rsb_mg_grammian, "PARRSB_RSB_MG_GRAMMIAN", 1);
  UPDATE_OPTION(rsb_mg_factor, "PARRSB_RSB_MG_FACTOR", 1);
  UPDATE_OPTION(rsb_mg_sagg, "PARRSB_RSB_MG_SMOOTH_AGGREGATION", 1);

#undef UPDATE_OPTION
}

static void print_options(const struct comm *c,
                          const parrsb_options *const options) {
#define PRINT_OPTION(OPT, STR, FMT)                                            \
  debug_print(c, options->verbose_level, "%s = " FMT "\n", STR, options->OPT)

  PRINT_OPTION(partitioner, "PARRSB_PARTITIONER", "%d");
  PRINT_OPTION(levels, "PARRSB_LEVELS", "%d");
  PRINT_OPTION(repair, "PARRSB_REPAIR", "%d");
  PRINT_OPTION(verbose_level, "PARRSB_VERBOSE_LEVEL", "%d");
  PRINT_OPTION(profile_level, "PARRSB_PROFILE_LEVEL", "%d");
  PRINT_OPTION(rsb_algo, "PARRSB_RSB_ALGO", "%d");
  PRINT_OPTION(rsb_pre, "PARRSB_RSB_PRE", "%d");
  PRINT_OPTION(rsb_max_iter, "PARRSB_RSB_MAX_ITER", "%d");
  PRINT_OPTION(rsb_max_passes, "PARRSB_RSB_MAX_PASSES", "%d");
  PRINT_OPTION(rsb_tol, "PARRSB_RSB_TOL", "%lf");
  PRINT_OPTION(rsb_dump_stats, "PARRSB_DUMP_STATS", "%d");
  PRINT_OPTION(rsb_mg_grammian, "PARRSB_RSB_MG_GRAMMIAN", "%d");
  PRINT_OPTION(rsb_mg_factor, "PARRSB_RSB_MG_FACTOR", "%d");
  PRINT_OPTION(rsb_mg_sagg, "PARRSB_RSB_MG_SMOOTH_AGGREGATION", "%d");

#undef PRINT_OPTION
}

static size_t load_balance(struct array *elist, uint nel, int nv,
                           const double *const coord,
                           const long long *const vtx, struct crystal *cr,
                           buffer *bfr) {
  struct comm *c = &cr->comm;
  slong out[2][1], wrk[2][1], in = nel;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0], nelg = out[1][0];

  uint nstar = nelg / c->np, nrem = nelg - nstar * c->np;
  slong lower = (nstar + 1) * nrem;

  size_t unit_size;
  if (vtx == NULL) // RCB
    unit_size = sizeof(struct rcb_element);
  else // RSB
    unit_size = sizeof(struct rsb_element);

  array_init_(elist, nel, unit_size, __FILE__, __LINE__);

  struct rcb_element *pe = (struct rcb_element *)calloc(1, unit_size);
  pe->origin = c->id;

  int ndim = (nv == 8) ? 3 : 2;
  for (uint e = 0; e < nel; ++e) {
    slong eg = pe->globalId = start + e + 1;
    if (nstar == 0)
      pe->proc = eg - 1;
    else if (eg <= lower)
      pe->proc = (eg - 1) / (nstar + 1);
    else
      pe->proc = (eg - 1 - lower) / nstar + nrem;

    pe->coord[0] = pe->coord[1] = pe->coord[2] = 0.0;
    for (int v = 0; v < nv; v++)
      for (int n = 0; n < ndim; n++)
        pe->coord[n] += coord[e * ndim * nv + v * ndim + n];
    for (int n = 0; n < ndim; n++)
      pe->coord[n] /= nv;

    array_cat_(unit_size, elist, pe, 1, __FILE__, __LINE__);
  }

  if (vtx != NULL) { // RSB
    struct rsb_element *pr = (struct rsb_element *)elist->ptr;
    for (uint e = 0; e < nel; e++) {
      for (int v = 0; v < nv; v++)
        pr[e].vertices[v] = vtx[e * nv + v];
    }
  }

  sarray_transfer_(elist, unit_size, offsetof(struct rcb_element, proc), 1, cr);
  if (vtx == NULL) // RCB
    sarray_sort(struct rcb_element, elist->ptr, elist->n, globalId, 1, bfr);
  else // RSB
    sarray_sort(struct rsb_element, elist->ptr, elist->n, globalId, 1, bfr);

  free(pe);
  return unit_size;
}

static void restore_original(int *part, int *seq, struct crystal *cr,
                             struct array *elist, size_t usize, buffer *bfr) {
  sarray_transfer_(elist, usize, offsetof(struct rcb_element, origin), 1, cr);
  uint nel = elist->n;

  if (usize == sizeof(struct rsb_element)) // RSB
    sarray_sort(struct rsb_element, elist->ptr, nel, globalId, 1, bfr);
  else if (usize == sizeof(struct rcb_element)) // RCB
    sarray_sort(struct rcb_element, elist->ptr, nel, globalId, 1, bfr);

  struct rcb_element *element;
  uint e;
  for (e = 0; e < nel; e++) {
    element = (struct rcb_element *)((char *)elist->ptr + e * usize);
    part[e] = element->origin; // element[e].origin;
  }

  if (seq != NULL) {
    for (e = 0; e < nel; e++) {
      element = (struct rcb_element *)((char *)elist->ptr + e * usize);
      seq[e] = element->seq; // element[e].seq;
    }
  }
}

int parrsb_part_mesh(int *part, int *seq, const long long *const vtx,
                     const double *const coord, const int nel, const int nv,
                     parrsb_options *const options, MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  update_options(options);

  const int verbose = options->verbose_level;
  {
    slong nelg = nel, wrk;
    comm_allreduce(&c, gs_long, gs_add, &nelg, 1, &wrk);
    debug_print(&c, verbose, "Running parRSB ..., nv = %d, nelg = %lld\n", nv,
                nelg);
  }

  print_options(&c, options);

  parrsb_barrier(&c);
  double t = comm_time();

  struct crystal cr;
  crystal_init(&cr, &c);

  buffer bfr;
  buffer_init(&bfr, (nel + 1) * sizeof(struct rsb_element));

  // Load balance input data
  debug_print(&c, verbose - 1, "Load balance: ...\n");
  struct array elist;
  size_t esize = load_balance(&elist, nel, nv, coord, vtx, &cr, &bfr);
  debug_print(&c, verbose - 1, "Load balance: done.\n");

  // Run RSB now
  debug_print(&c, verbose - 1, "Running partitioner: ...\n");
  struct comm ca;
  comm_split(&c, elist.n > 0, c.id, &ca);
  metric_init();
  if (elist.n > 0) {
    slong out[2][1], wrk[2][1], in = elist.n;
    comm_scan(out, &ca, gs_long, gs_add, &in, 1, wrk);
    slong nelg = out[1][0];

    int ndim = (nv == 8) ? 3 : 2;
    switch (options->partitioner) {
    case 0:
      rsb(&elist, nv, options, &ca, &bfr);
      break;
    case 1:
      rcb(&elist, esize, ndim, &ca, &bfr);
      break;
    case 2:
      rib(&elist, esize, ndim, &ca, &bfr);
      break;
    default:
      break;
    }

    metric_rsb_print(&ca, options->profile_level);
  }
  metric_finalize(), comm_free(&ca);
  debug_print(&c, verbose - 1, "Running partitioner: done.\n");

  debug_print(&c, verbose - 1, "Restore original input: ...\n");
  restore_original(part, seq, &cr, &elist, esize, &bfr);
  debug_print(&c, verbose - 1, "Restore original input: done.\n");

  // Report time and finish
  debug_print(&c, verbose, "par%s finished in %g seconds.\n",
              ALGO[options->partitioner], comm_time() - t);

  array_free(&elist), buffer_free(&bfr), crystal_free(&cr), comm_free(&c);

  return 0;
}

void fparrsb_part_mesh(int *part, int *seq, long long *vtx, double *coord,
                       int *nel, int *nv, int *options, int *comm, int *err) {
  *err = 1;
  MPI_Comm c = MPI_Comm_f2c(*comm);
  parrsb_options opt = parrsb_default_options;
  *err = parrsb_part_mesh(part, seq, vtx, coord, *nel, *nv, &opt, c);
}

int parrsb_part_mesh_v2(int *part, const long long *const vtx,
                        const double *const coord, const int nel, const int nv,
                        const int *const tag, parrsb_options *const options,
                        MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  update_options(options);

  const int verbose = options->verbose_level;
  {
    slong nelg = nel, wrk;
    comm_allreduce(&c, gs_long, gs_add, &nelg, 1, &wrk);
    const int verbose = options->verbose_level;
    debug_print(&c, verbose, "Running parRSB v2..., nv = %d, nelg = %lld\n", nv,
                nelg);
  }

  print_options(&c, options);

  struct tag_t {
    uint p, tag, seq, tagn;
  };

  struct array tags;
  array_init(struct tag_t, &tags, nel);

  buffer bfr;
  buffer_init(&bfr, nel * sizeof(struct tag_t));

  {
    struct tag_t tt;
    for (uint i = 0; i < nel; i++) {
      tt.seq = i, tt.tag = tag[i], tt.p = tt.tag % c.np;
      array_cat(struct tag_t, &tags, &tt, 1);
    }
    sarray_sort(struct tag_t, tags.ptr, tags.n, tag, 0, &bfr);
  }

  struct array unique;
  array_init(struct tag_t, &unique, 1024);

  if (tags.n > 0) {
    const struct tag_t *const pt = (const struct tag_t *const)tags.ptr;
    array_cat(struct tag_t, &unique, &pt[0], 1);
    for (uint i = 1; i < tags.n; i++) {
      if (pt[i].tag > pt[i - 1].tag)
        array_cat(struct tag_t, &unique, &pt[i], 1);
    }
  }

  struct crystal cr;
  crystal_init(&cr, &c);

  sint out[2][1];
  {
    sarray_transfer(struct tag_t, &unique, p, 1, &cr);
    sarray_sort(struct tag_t, unique.ptr, unique.n, tag, 0, &bfr);

    const struct tag_t *const pu = (const struct tag_t *const)unique.ptr;
    sint in = 0;
    if (unique.n > 0) {
      in = 1;
      for (uint i = 1; i < unique.n; i++) {
        if (pu[i].tag > pu[i - 1].tag)
          in++;
      }
    }

    sint wrk[2][1];
    comm_scan(out, &c, gs_int, gs_add, &in, 1, wrk);
  }
  const uint num_tags = out[1][0], tag_start = out[0][0];

  debug_print(&c, verbose, "\tNum tags: %d\n", num_tags);
  if (c.np % num_tags != 0) {
    if (c.id == 0) {
      fprintf(stderr,
              "Number of processes must be a multiple of number of tags: "
              "processes = %d, tags = %d.\n",
              c.np, num_tags);
    }
    exit(EXIT_FAILURE);
  }

  {
    struct tag_t *const pu = (struct tag_t *const)unique.ptr;
    uint start = tag_start;
    if (unique.n > 0) {
      pu[0].tagn = start;
      for (uint i = 1; i < unique.n; i++) {
        if (pu[i].tag > pu[i - 1].tag)
          start++;
        pu[i].tagn = start;
      }
    }

    sarray_transfer(struct tag_t, &unique, p, 0, &cr);
    sarray_sort(struct tag_t, unique.ptr, unique.n, tag, 0, &bfr);
  }

  const uint chunk_size = c.np / num_tags;
  debug_print(&c, verbose, "\tProcesses per tag: %d\n", chunk_size);
  {
    struct tag_t *const pt = (struct tag_t *const)tags.ptr;
    const struct tag_t *const pu = (const struct tag_t *const)unique.ptr;
    for (uint i = 0, s = 0; i < unique.n; i++) {
      uint e = s + 1;
      assert(pt[s].tag == pu[i].tag);
      while (e < tags.n && pt[e].tag == pu[i].tag)
        e++;
      for (uint j = s; j < e; j++)
        pt[j].p = chunk_size * pu[i].tagn + pt[i].seq % chunk_size;
      s = e;
    }

    sarray_sort(struct tag_t, tags.ptr, tags.n, seq, 0, &bfr);
  }

  struct element_t {
    uint proc, part, seq;
    scalar coord[MAXDIM * MAXNV];
    slong vertices[MAXNV];
  };

  struct array elements;
  array_init(struct element_t, &elements, nel);

  debug_print(&c, verbose - 1,
              "\tPack element data for transfer. tags.n=%u, nel=%u\n", tags.n,
              nel);
  const int ndim = (nv == 8) ? 3 : 2;
  {
    assert(tags.n == nel);
    const struct tag_t *const pt = (const struct tag_t *const)tags.ptr;
    struct element_t et;
    for (uint i = 0; i < tags.n; i++) {
      et.proc = pt[i].p, et.seq = i;
      for (uint j = 0; j < nv; j++) {
        et.vertices[j] = vtx[i * nv + j];
        for (uint k = 0; k < ndim; k++)
          et.coord[j * ndim + k] = coord[i * nv * ndim + j * ndim + k];
      }
      array_cat(struct element_t, &elements, &et, 1);
    }

    sarray_transfer(struct element_t, &elements, proc, 1, &cr);
  }

  debug_print(&c, verbose - 1, "\tCopy element data for feeding to parRSB.\n");
  long long *lvtx = tcalloc(long long, (elements.n + 1) * nv);
  double *lcoord = tcalloc(double, (elements.n + 1) * nv * ndim);
  {
    const struct element_t *const pe =
        (const struct element_t *const)elements.ptr;
    for (uint e = 0; e < elements.n; e++) {
      for (uint j = 0; j < nv; j++) {
        lvtx[e * nv + j] = pe[e].vertices[j];
        for (uint k = 0; k < ndim; k++)
          lcoord[e * nv * ndim + j * ndim + k] = pe[e].coord[j * ndim + k];
      }
    }
  }

  debug_print(&c, verbose - 1, "\tRun parRSB locally in a tag now.\n");
  {
    int *lpart = tcalloc(int, elements.n + 1);

    MPI_Comm local;
    MPI_Comm_split(comm, c.id / chunk_size, c.id, &local);
    options->verbose_level = 0;
    options->profile_level = 0;
    parrsb_part_mesh(lpart, NULL, lvtx, lcoord, elements.n, nv, options, local);
    MPI_Comm_free(&local);

    struct element_t *const pe = (struct element_t *const)elements.ptr;
    for (uint e = 0; e < elements.n; e++) {
      pe[e].part = lpart[e] + (c.id / chunk_size) * chunk_size;
      assert(pe[e].part >= 0 && pe[e].part < c.np);
    }
    free(lpart);

    sarray_transfer(struct element_t, &elements, proc, 0, &cr);
    assert(nel == elements.n);
  }
  free(lvtx), free(lcoord);

  {
    sarray_sort(struct element_t, elements.ptr, elements.n, seq, 0, &bfr);
    const struct element_t *const pe =
        (const struct element_t *const)elements.ptr;
    for (uint i = 0; i < nel; i++)
      part[i] = pe[i].part;
  }

  array_free(&elements), array_free(&unique), array_free(&tags);
  buffer_free(&bfr), crystal_free(&cr), comm_free(&c);

  return 0;
}

void fparrsb_partmesh_v2(int *part, const long long *const vtx,
                         const double *const coord, const int *const nel,
                         const int *const nv, const int *const tag,
                         const int *const options, const int *const comm,
                         int *err) {
  *err = 1;
  MPI_Comm c = MPI_Comm_f2c(*comm);
  parrsb_options opt = parrsb_default_options;
  *err = parrsb_part_mesh_v2(part, vtx, coord, *nel, *nv, tag, &opt, c);
}
