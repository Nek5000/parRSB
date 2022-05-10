#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <genmap-impl.h>
#include <parRSB.h>

parrsb_options parrsb_default_options = {0, 2, 1, 0, 1, 0, 4, 1};

static char *ALGO[3] = {"RSB", "RCB", "RIB"};

static int if_number(const char *c) {
  int i;
  for (i = 0; i < strnlen(c, 32) && isdigit(c[i]); i++)
    ;
  return i == strnlen(c, 32);
}

#define UPDATE_OPTION(opt, str)                                                \
  {                                                                            \
    const char *val = getenv(str);                                             \
    if (val != NULL && if_number(val))                                         \
      options->opt = atoi(val);                                                \
  }

static void update_options(parrsb_options *options) {
  UPDATE_OPTION(partitioner, "PARRSB_PARTITIONER");
  UPDATE_OPTION(verbose_level, "PARRSB_VERBOSE_LEVEL");
  UPDATE_OPTION(profile_level, "PARRSB_PROFILE_LEVEL");
  UPDATE_OPTION(rsb_algo, "PARRSB_RSB_ALGO");
  UPDATE_OPTION(rsb_pre, "PARRSB_RSB_PRE");
  UPDATE_OPTION(rsb_grammian, "PARRSB_RSB_GRAMMIAN");
  UPDATE_OPTION(repair, "PARRSB_REPAIR");
  UPDATE_OPTION(rsb_mg_factor, "PARRSB_RSB_MG_FACTOR");
}

#undef UPDATE_OPTION
#define PRINT_OPTION(opt, str) printf("%s = %d\n", str, options->opt)

static void print_options(parrsb_options *options) {
  PRINT_OPTION(partitioner, "PARRSB_PARTITIONER");
  PRINT_OPTION(verbose_level, "PARRSB_VERBOSE_LEVEL");
  PRINT_OPTION(profile_level, "PARRSB_PROFILE_LEVEL");
  PRINT_OPTION(rsb_algo, "PARRSB_RSB_ALGO");
  PRINT_OPTION(rsb_pre, "PARRSB_RSB_PRE");
  PRINT_OPTION(rsb_grammian, "PARRSB_RSB_GRAMMIAN");
  PRINT_OPTION(repair, "PARRSB_REPAIR");
  PRINT_OPTION(rsb_mg_factor, "PARRSB_RSB_MG_FACTOR");
}

#undef PRINT_OPTION

// Load balance input data
static size_t load_balance(struct array *elist, uint nel, int nv, double *coord,
                           long long *vtx, struct crystal *cr, buffer *bfr) {
  slong out[2][1], buf[2][1];
  slong in = nel;
  comm_scan(out, &cr->comm, gs_long, gs_add, &in, 1, buf);
  slong start = out[0][0];
  slong nelg = out[1][0];

  int size = cr->comm.np;
  uint nstar = nelg / size;
  uint nrem = nelg - nstar * size;
  slong lower = (nstar + 1) * nrem;

  size_t unit_size;
  struct rcb_element *element = NULL;

  if (vtx == NULL) // RCB
    unit_size = sizeof(struct rcb_element);
  else
    unit_size = sizeof(struct rsb_element);

  element = (struct rcb_element *)calloc(1, unit_size);
  element->origin = cr->comm.id;

  array_init_(elist, nel, unit_size, __FILE__, __LINE__);

  int ndim = (nv == 8) ? 3 : 2;
  int e, n, v;
  for (e = 0; e < nel; ++e) {
    slong eg = element->globalId = start + e + 1;
    if (eg <= lower)
      element->proc = (eg - 1) / (nstar + 1);
    else if (nstar != 0)
      element->proc = (eg - 1 - lower) / nstar + nrem;

    element->coord[0] = element->coord[1] = element->coord[2] = 0.0;
    for (v = 0; v < nv; v++)
      for (n = 0; n < ndim; n++)
        element->coord[n] += coord[e * ndim * nv + v * ndim + n];
    for (n = 0; n < ndim; n++)
      element->coord[n] /= nv;

    array_cat_(unit_size, elist, element, 1, __FILE__, __LINE__);
  }
  // Sanity check
  assert(elist->n == nel);

  if (vtx != NULL) { // RSB
    struct rsb_element *elements = elist->ptr;
    for (e = 0; e < nel; e++)
      for (v = 0; v < nv; v++)
        elements[e].vertices[v] = vtx[e * nv + v];
  }

  sarray_transfer_(elist, unit_size, offsetof(struct rcb_element, proc), 1, cr);
  nel = elist->n;
  if (vtx != NULL) // RSB
    sarray_sort(struct rsb_element, elist->ptr, nel, globalId, 1, bfr);
  else
    sarray_sort(struct rcb_element, elist->ptr, nel, globalId, 1, bfr);

  free(element);

  return unit_size;
}

static void restore_original(int *part, int *seq, struct crystal *cr,
                             struct array *elist, size_t unit_size,
                             buffer *bfr) {
  sarray_transfer_(elist, unit_size, offsetof(struct rcb_element, origin), 1,
                   cr);
  uint nel = elist->n;

  if (unit_size == sizeof(struct rsb_element)) // RSB
    sarray_sort(struct rsb_element, elist->ptr, nel, globalId, 1, bfr);
  else if (unit_size == sizeof(struct rcb_element)) // RCB
    sarray_sort(struct rcb_element, elist->ptr, nel, globalId, 1, bfr);

  struct rcb_element *element;
  uint e;
  for (e = 0; e < nel; e++) {
    element = (struct rcb_element *)((char *)elist->ptr + e * unit_size);
    part[e] = element->origin; // element[e].origin;
  }

  if (seq != NULL) {
    for (e = 0; e < nel; e++) {
      element = (struct rcb_element *)((char *)elist->ptr + e * unit_size);
      seq[e] = element->seq; // element[e].seq;
    }
  }
}

/*
 * part = [nel], out,
 * seq = [nel], out,
 * vtx = [nel x nv], in,
 * coord = [nel x nv x ndim], in,
 * nel = in,
 * nv = in,
 * options = in */
int parrsb_part_mesh(int *part, int *seq, long long *vtx, double *coord,
                     int nel, int nv, parrsb_options options, MPI_Comm comm) {
  update_options(&options);

  struct comm c;
  comm_init(&c, comm);
  if (c.id == 0) {
    if (options.verbose_level > 0)
      printf("Running parRSB ...\n");
    if (options.verbose_level > 1)
      print_options(&options);
    fflush(stdout);
  }

  genmap_barrier(&c);
  double t = comm_time();

#if defined(GENMAP_OCCA)
  occa_init("CUDA", 0, 0, &c);
#endif

  struct crystal cr;
  crystal_init(&cr, &c);

  buffer bfr;
  buffer_init(&bfr, nel * sizeof(struct rsb_element));

  // Load balance input data
  struct array elist;
  size_t elem_size = load_balance(&elist, nel, nv, coord, vtx, &cr, &bfr);

  // Run RSB now
  MPI_Comm comma;
  MPI_Comm_split(c.c, nel > 0, c.id, &comma);

  if (nel > 0) {
    metric_init();

    struct comm ca;
    comm_init(&ca, comma);

    slong out[2][1], buf[2][1];
    slong in = nel;
    comm_scan(out, &ca, gs_long, gs_add, &in, 1, buf);
    slong nelg = out[1][0];

    if (ca.np > nelg) {
      if (ca.id == 0)
        printf("Total number of elements is smaller than the "
               "number of processors.\n"
               "Run with smaller number of processors.\n");
      return 1;
    }

    int ndim = (nv == 8) ? 3 : 2;
    switch (options.partitioner) {
    case 0:
      rsb(&elist, &options, nv, &ca, &bfr);
      break;
    case 1:
      rcb(&elist, elem_size, ndim, &ca, &bfr);
      break;
    case 2:
      rib(&elist, elem_size, ndim, &ca, &bfr);
      break;
    default:
      break;
    }

    metric_print(&ca, options.profile_level);

    comm_free(&ca);
    metric_finalize();
  }

#if defined(GENMAP_OCCA)
  occa_free();
#endif

  MPI_Comm_free(&comma);

  restore_original(part, seq, &cr, &elist, elem_size, &bfr);

  // Report time and finish
  genmap_barrier(&c);
  t = comm_time() - t;
  if (c.id == 0 && options.verbose_level > 0) {
    printf("par%s finished in %g s\n", ALGO[options.partitioner], t);
    fflush(stdout);
  }

  array_free(&elist);
  buffer_free(&bfr);
  crystal_free(&cr);
  comm_free(&c);

  return 0;
}

void fparrsb_part_mesh(int *part, int *seq, long long *vtx, double *coord,
                       int *nel, int *nv, int *options, int *comm, int *err) {
  *err = 1;
  comm_ext c = MPI_Comm_f2c(*comm);
  // TODO: Convert int options to parrsb_options instead of default options
  parrsb_options opt = parrsb_default_options;
  *err = parrsb_part_mesh(part, seq, vtx, coord, *nel, *nv, opt, c);
}
