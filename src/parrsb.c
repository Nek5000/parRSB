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
                           const double *const xyz, const long long *const vtx,
                           struct crystal *cr, buffer *bfr) {
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
        pe->coord[n] += xyz[e * ndim * nv + v * ndim + n];
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

static void restore_original(int *part, struct crystal *cr, struct array *elist,
                             size_t usize, buffer *bfr) {
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
}

void parrsb_part_mesh_v1(int *part, const long long *const vtx,
                         const double *const xyz, const int nel, const int nv,
                         const parrsb_options *const options,
                         const struct comm *const c, struct crystal *const cr,
                         buffer *const bfr) {
  const int verbose = options->verbose_level - 1;
  debug_print(c, verbose, "Load balance ...\n");

  struct array elist;
  size_t esize = load_balance(&elist, nel, nv, xyz, vtx, cr, bfr);

  debug_print(c, verbose, "Running partitioner ...\n");
  struct comm ca;
  comm_split(c, elist.n > 0, c->id, &ca);
  if (elist.n > 0) {
    slong out[2][1], wrk[2][1], in = elist.n;
    comm_scan(out, &ca, gs_long, gs_add, &in, 1, wrk);
    slong nelg = out[1][0];

    int ndim = (nv == 8) ? 3 : 2;
    switch (options->partitioner) {
    case 0:
      rsb(&elist, nv, options, &ca, bfr);
      break;
    case 1:
      rcb(&elist, esize, ndim, &ca, bfr);
      break;
    case 2:
      rib(&elist, esize, ndim, &ca, bfr);
      break;
    default:
      break;
    }
  }
  comm_free(&ca);

  debug_print(c, verbose, "Restore original input: ...\n");
  restore_original(part, cr, &elist, esize, bfr);

  array_free(&elist);
}

void parrsb_check_tagged_partitions(const long long *const eids,
                                    const long long *const vtx, const int nel,
                                    const int nv, const int ntags,
                                    const struct comm *const c) {
  // Check if the input elements are sorted by global id.
  {
    sint sorted = 1;
    for (uint i = 1; i < nel; i++) {
      if (eids[i] < eids[i - 1]) {
        sorted = 0;
        break;
      }
    }

    sint wrk;
    comm_allreduce(c, gs_int, gs_min, &sorted, 1, &wrk);
    if (!sorted) {
      if (c->id == 0) {
        fprintf(stderr, "Input elements are not sorted.\n");
        fflush(stderr);
      }
      exit(EXIT_FAILURE);
    }
  }

  // Number the elements within the each tag id and setup a gs handle based on
  // 2D element id.
  const uint tag_id = c->id / ntags;
  struct comm lc;
  struct gs_data *gse = NULL;
  {
    comm_split(c, tag_id, c->id, &lc);

    slong out[2][1], wrk[2][1], in = nel;
    comm_scan(out, &lc, gs_long, gs_add, &in, 1, wrk);
    slong start = out[0][0];

    slong *lids = tcalloc(slong, nel);
    for (uint i = 0; i < nel; i++)
      lids[i] = start + i;

    gse = gs_setup(lids, nel, c, 0, gs_pairwise, 0);
    free(lids);
  }

  // Setup a local gs handle based on the original gs vertex ids.
  const size_t size = nel * nv;
  buffer bfr;
  buffer_init(&bfr, size);
  sint *mul = tcalloc(sint, size);
  {
    struct gs_data *gsl = gs_setup(vtx, size, &lc, 0, gs_pairwise, 0);
    for (uint i = 0; i < size; i++)
      mul[i] = 1;
    gs(mul, gs_int, gs_add, 0, gsl, &bfr);
    gs_free(gsl);
  }

  // Now let's compare the multiplicity across the layers.
  {
    sint *lmin = tcalloc(sint, nel);
    sint *lmax = tcalloc(sint, nel);
    for (uint v = 0; v < nv; v++) {
      for (uint e = 0; e < nel; e++) {
        lmin[e] = mul[e * nv + v];
        lmax[e] = mul[e * nv + v];
      }

      gs(lmin, gs_int, gs_min, 0, gse, &bfr);
      gs(lmax, gs_int, gs_max, 0, gse, &bfr);

      for (uint e = 0; e < nel; e++)
        assert(lmin[e] == lmax[e]);
    }

    free(lmin), free(lmax);
  }

  free(mul);
  buffer_free(&bfr);
  gs_free(gse);
  comm_free(&lc);

  return;
}

void parrsb_part_mesh_v2(int *part, const long long *const vtx,
                         const double *const xyz, const int *const tag,
                         const int nel, const int nv,
                         parrsb_options *const options,
                         const struct comm *const c, struct crystal *const cr,
                         buffer *const bfr) {
  const int verbose = options->verbose_level - 1;
  debug_print(c, verbose, "Find number of tags in the mesh ...\n");

  struct tag_t {
    uint p, tag, seq, tagn;
  };

  struct array tags;
  array_init(struct tag_t, &tags, nel);

  {
    struct tag_t tt;
    for (uint i = 0; i < nel; i++) {
      tt.seq = i, tt.tag = tag[i], tt.p = tt.tag % c->np;
      array_cat(struct tag_t, &tags, &tt, 1);
    }
    sarray_sort(struct tag_t, tags.ptr, tags.n, tag, 0, bfr);
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

  sint out[2][1];
  {
    sarray_transfer(struct tag_t, &unique, p, 1, cr);
    sarray_sort(struct tag_t, unique.ptr, unique.n, tag, 0, bfr);

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
    comm_scan(out, c, gs_int, gs_add, &in, 1, wrk);
  }
  const uint num_tags = out[1][0], tag_start = out[0][0];

  debug_print(c, verbose, "Num tags: %d\n", num_tags);
  if (c->np % num_tags != 0) {
    if (c->id == 0) {
      fprintf(stderr,
              "Number of processes must be a multiple of number of tags: "
              "processes = %d, tags = %d.\n",
              c->np, num_tags);
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

    sarray_transfer(struct tag_t, &unique, p, 0, cr);
    sarray_sort(struct tag_t, unique.ptr, unique.n, tag, 0, bfr);
  }

  const uint chunk_size = c->np / num_tags;
  debug_print(c, verbose, "Processes per tag: %d\n", chunk_size);
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

    sarray_sort(struct tag_t, tags.ptr, tags.n, seq, 0, bfr);
  }
  array_free(&unique);

  struct element_t {
    uint proc, part, seq;
    scalar xyz[MAXDIM * MAXNV];
    slong vertices[MAXNV];
  };

  struct array elements;
  array_init(struct element_t, &elements, nel);

  debug_print(c, verbose,
              "Pack element data for transfering. tags.n=%u, nel=%u\n", tags.n,
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
          et.xyz[j * ndim + k] = xyz[i * nv * ndim + j * ndim + k];
      }
      array_cat(struct element_t, &elements, &et, 1);
    }

    sarray_transfer(struct element_t, &elements, proc, 1, cr);
  }
  array_free(&tags);

  debug_print(c, verbose, "Copy element data for feeding to parRSB.\n");
  long long *lvtx = tcalloc(long long, (elements.n + 1) * nv);
  double *lxyz = tcalloc(double, (elements.n + 1) * nv * ndim);
  {
    const struct element_t *const pe =
        (const struct element_t *const)elements.ptr;
    for (uint e = 0; e < elements.n; e++) {
      for (uint j = 0; j < nv; j++) {
        lvtx[e * nv + j] = pe[e].vertices[j];
        for (uint k = 0; k < ndim; k++)
          lxyz[e * nv * ndim + j * ndim + k] = pe[e].xyz[j * ndim + k];
      }
    }
  }

  debug_print(c, verbose, "Run parRSB locally within a tag now.\n");
  {
    int *lpart = tcalloc(int, elements.n + 1);

    struct comm lc;
    comm_split(c, c->id / chunk_size, c->id, &lc);

    struct crystal lcr;
    crystal_init(&lcr, &lc);

    options->verbose_level = 0;
    options->profile_level = 0;
    parrsb_part_mesh_v1(lpart, lvtx, lxyz, elements.n, nv, options, &lc, &lcr,
                        bfr);
    crystal_free(&lcr), comm_free(&lc);

    struct element_t *const pe = (struct element_t *const)elements.ptr;
    for (uint e = 0; e < elements.n; e++) {
      pe[e].part = lpart[e] + (c->id / chunk_size) * chunk_size;
      assert(pe[e].part >= 0 && pe[e].part < c->np);
    }
    free(lpart);

    sarray_transfer(struct element_t, &elements, proc, 0, cr);
    assert(nel == elements.n);
  }
  free(lvtx), free(lxyz);

  {
    sarray_sort(struct element_t, elements.ptr, elements.n, seq, 0, bfr);
    const struct element_t *const pe =
        (const struct element_t *const)elements.ptr;
    for (uint i = 0; i < nel; i++)
      part[i] = pe[i].part;
  }

  array_free(&elements);
}

int parrsb_part_mesh(int *part, const long long *const vtx,
                     const double *const xyz, const int *const tag,
                     const int nel, const int nv, parrsb_options *const options,
                     MPI_Comm comm) {
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

  buffer bfr;
  buffer_init(&bfr, (nel + 1) * 72);

  struct crystal cr;
  crystal_init(&cr, &c);

  parrsb_barrier(&c);
  const double t = comm_time();

  metric_init();
  const uint profile_level = options->profile_level;
  if (tag != NULL)
    parrsb_part_mesh_v2(part, vtx, xyz, tag, nel, nv, options, &c, &cr, &bfr);
  else
    parrsb_part_mesh_v1(part, vtx, xyz, nel, nv, options, &c, &cr, &bfr);
  metric_rsb_print(&c, profile_level);
  metric_finalize();

  debug_print(&c, verbose, "par%s finished in %g seconds.\n",
              ALGO[options->partitioner], comm_time() - t);

  crystal_free(&cr);
  buffer_free(&bfr);
  comm_free(&c);
}
