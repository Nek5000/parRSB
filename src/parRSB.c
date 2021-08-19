#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <genmap-impl.h>
#include <parRSB.h>

parrsb_options parrsb_default_options = {0, 0, 0, 0, 1, 1, 1};

static int if_number(const char *c) {
  int i;
  for (i = 0; i < strlen(c); i++)
    if (!isdigit(c[i]))
      return 0;
  return 1;
}

#define init_option(opt, str)                                                  \
  do {                                                                         \
    const char *val = getenv(str);                                             \
    if (val != NULL && if_number(val))                                         \
      options->opt = atoi(val);                                                \
  } while (0)

static void init_options(parrsb_options *options) {
  init_option(partitioner, "PARRSB_PARTITIONER");
  init_option(debug_level, "PARRSB_DEBUG_LEVEL");
  init_option(profile_level, "PARRSB_PROFILE_LEVEL");
  init_option(rsb_algo, "PARRSB_RSB_ALGO");
  init_option(rsb_pre, "PARRSB_RSB_PRE");
  init_option(rsb_grammian, "PARRSB_RSB_GRAMMIAN");
  init_option(repair, "PARRSB_REPAIR");
}

#undef init_option

#define print_option(opt, str) printf("%s = %d\n", str, options->opt)

static void print_options(parrsb_options *options) {
  print_option(partitioner, "PARRSB_PARTITIONER");
  print_option(debug_level, "PARRSB_DEBUG_LEVEL");
  print_option(profile_level, "PARRSB_PROFILE_LEVEL");
  print_option(rsb_algo, "PARRSB_RSB_ALGO");
  print_option(rsb_pre, "PARRSB_RSB_PRE");
  print_option(rsb_grammian, "PARRSB_RSB_GRAMMIAN");
  print_option(repair, "PARRSB_REPAIR");
}

#undef print_option

/*
 * part = [nel], out,
 * seq = [nel], out,
 * vtx = [nel x nv], in,
 * coord = [nel x nv x ndim], in,
 * nel = in,
 * nv = in,
 * options = in */
int parRSB_partMesh(int *part, int *seq, long long *vtx, double *coord, int nel,
                    int nv, parrsb_options options, MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);
  int rank = c.id;
  int size = c.np;

  if (rank == 0) {
    printf("Running parRSB ...\n");
    fflush(stdout);
  }

  comm_barrier(&c);
  double t = comm_time();

  init_options(&options);
  if (rank == 0) {
    print_options(&options);
    fflush(stdout);
  }

  struct crystal cr;
  crystal_init(&cr, &c);

  buffer bfr;
  buffer_init(&bfr, 1024);

  /* Load balance input data */
  struct array elist;
  genmap_load_balance(&elist, nel, nv, coord, vtx, &cr, &bfr);

  /* Run RSB now */
  comm_ext comm_rsb;
#ifdef MPI
  MPI_Comm_split(c.c, nel > 0, rank, &comm_rsb);
#endif

  // TODO: Move this into another file
  if (nel > 0) {
    metric_init();

    genmap_handle h;
    genmap_init(&h, comm_rsb, &options);

    genmap_set_elements(h, &elist);
    genmap_comm_scan(h, genmap_global_comm(h));
    genmap_set_nvertices(h, nv);

    GenmapLong nelg = genmap_get_partition_nel(h);
    GenmapInt id = genmap_comm_rank(genmap_global_comm(h));
    GenmapInt size_ = genmap_comm_size(genmap_global_comm(h));

    if (size_ > nelg) {
      if (id == 0)
        printf("Total number of elements is smaller than the "
               "number of processors.\n"
               "Run with smaller number of processors.\n");
      return 1;
    }

    switch (options.partitioner) {
    case 0:
      genmap_rsb(h);
      break;
    case 1:
      genmap_rcb(h);
      break;
    case 2:
      genmap_rib(h);
      break;
    default:
      break;
    }

    genmap_finalize(h);

    if (options.debug_level > 0)
      metric_print(&c);
    metric_finalize();
  }

#ifdef MPI
  MPI_Comm_free(&comm_rsb);
#endif

  genmap_restore_original(part, seq, &cr, &elist, &bfr);

  /* Report time and finish */
  comm_barrier(&c);
  t = comm_time() - t;
  if (rank == 0) {
    printf("parRSB finished in %g s\n", t);
    fflush(stdout);
  }

  array_free(&elist);
  buffer_free(&bfr);
  crystal_free(&cr);
  comm_free(&c);

  return 0;
}

void fparRSB_partMesh(int *part, int *seq, long long *vtx, double *coord,
                      int *nel, int *nv, int *options, int *comm, int *err) {
  *err = 1;
  comm_ext c = MPI_Comm_f2c(*comm);
  /* TODO: Convert int options to parrsb_options instead of default options */
  parrsb_options opt = parrsb_default_options;
  *err = parRSB_partMesh(part, seq, vtx, coord, *nel, *nv, opt, c);
}
