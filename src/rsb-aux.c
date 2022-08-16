#include "coarse.h"
#include "genmap-impl.h"
#include "ilu.h"
#include "metrics.h"
#include "sort.h"
#include <limits.h>
#include <time.h>

extern int rcb(struct array *elements, size_t unit_size, int ndim,
               struct comm *c, buffer *bfr);
extern int fiedler(struct array *elements, int nv, parrsb_options *options,
                   struct comm *gsc, buffer *buf);
extern int repair_partitions(struct array *elements, int nv, struct comm *tc,
                             struct comm *lc, int bin, struct comm *gc,
                             buffer *bfr);
extern int balance_partitions(struct array *elements, int nv, struct comm *lc,
                              struct comm *gc, int bin, buffer *bfr);

static void check_rsb_partition(struct comm *gc, parrsb_options *opts) {
  int max_levels = log2ll(gc->np);
  int miter = opts->rsb_max_iter, mpass = opts->rsb_lanczos_max_restarts;

  for (int i = 0; i < max_levels; i++) {
    sint converged = 1;
    int val = (int)metric_get_value(i, RSB_FIEDLER_INNER_NITER);
    if (opts->rsb_algo == 0) {
      if (val >= miter * mpass)
        converged = 0;
    } else if (opts->rsb_algo == 1) {
      if (val >= miter)
        converged = 0;
    }

    sint ibfr;
    comm_allreduce(gc, gs_int, gs_min, &converged, 1, &ibfr);

    if (converged == 0) {
      double final = (double)metric_get_value(i, TOL_FNL), dbfr;
      comm_allreduce(gc, gs_double, gs_min, &final, 1, &dbfr);

      double target = (double)metric_get_value(i, TOL_TGT);
      comm_allreduce(gc, gs_double, gs_min, &target, 1, &dbfr);

      if (gc->id == 0) {
        printf("Warning: Partitioner only reached a tolerance of %lf given %lf "
               "after %d x %d iterations in Level=%d!\n",
               final, target, mpass, miter, i);
        fflush(stdout);
      }
    }

    sint ncomps, minc, maxc;
    ncomps = minc = maxc = (sint)metric_get_value(i, RSB_COMPONENTS);
    comm_allreduce(gc, gs_int, gs_min, &minc, 1, &ibfr);
    comm_allreduce(gc, gs_int, gs_max, &maxc, 1, &ibfr);

    if (maxc > 1 && gc->id == 0) {
      printf("Warning: Partition created %d/%d (min/max) disconnected "
             "components in Level=%d!\n",
             minc, maxc, i);
      fflush(stdout);
    }
  }
}

static void get_part(sint *np, sint *nid, int two_lvl, struct comm *lc,
                     struct comm *nc) {
  if (two_lvl) {
    sint out[2][1], wrk[2][1], in = (nc->id == 0);
    comm_scan(out, lc, gs_int, gs_add, &in, 1, &wrk);
    *nid = (nc->id == 0) * out[0][0], *np = out[1][0];
    comm_allreduce(nc, gs_int, gs_max, nid, 1, wrk);
  } else {
    *np = lc->np, *nid = lc->id;
  }
}

int rsb(struct array *elements, int nv, parrsb_options *options,
        struct comm *gc, buffer *bfr) {
  int ndim = (nv == 8) ? 3 : 2;

  struct comm lc, tc;
  comm_dup(&lc, gc);

  struct comm nc;
  if (options->two_level) {
#ifdef MPI
    MPI_Comm node;
    MPI_Comm_split_type(lc.c, MPI_COMM_TYPE_SHARED, lc.id, MPI_INFO_NULL,
                        &node);
    comm_init(&nc, node);
    MPI_Comm_free(&node);
#else
    comm_init(&nc, 1);
#endif
  }

  sint np, nid;
  get_part(&np, &nid, options->two_level, &lc, &nc);
  if (options->two_level && options->verbose_level) {
    if (gc->id == 0)
      printf("Number of nodes = %d\n", np);
    fflush(stdout);
  }

  while (np > 1) {
    // Run RCB, RIB pre-step or just sort by global id
    metric_tic(&lc, RSB_PRE);
    switch (options->rsb_pre) {
    case 0: // Sort by global id
      parallel_sort(struct rsb_element, elements, globalId, gs_long, 0, 1, &lc,
                    bfr);
      break;
    case 1: // RCB
      rcb(elements, sizeof(struct rsb_element), ndim, &lc, bfr);
      break;
    case 2: // RIB
      rib(elements, sizeof(struct rsb_element), ndim, &lc, bfr);
      break;
    default:
      break;
    }
    metric_toc(&lc, RSB_PRE);

    // Run fiedler
    metric_tic(&lc, RSB_FIEDLER);
    fiedler(elements, nv, options, &lc, bfr);
    metric_toc(&lc, RSB_FIEDLER);

    // Sort by Fiedler vector
    metric_tic(&lc, RSB_SORT);
    parallel_sort_2(struct rsb_element, elements, fiedler, gs_double, globalId,
                    gs_long, 0, 1, &lc, bfr);
    metric_toc(&lc, RSB_SORT);

    // Bisect, repair and balance
    metric_tic(&lc, RSB_REPAIR_BALANCE);
    int bin = (nid >= (np + 1) / 2);
    comm_split(&lc, bin, lc.id, &tc);
    if (options->repair > 0)
      repair_partitions(elements, nv, &tc, &lc, bin, gc, bfr);
    balance_partitions(elements, nv, &tc, &lc, bin, bfr);
    metric_toc(&lc, RSB_REPAIR_BALANCE);

    comm_free(&lc);
    comm_dup(&lc, &tc);
    comm_free(&tc);

    get_part(&np, &nid, options->two_level, &lc, &nc);
    metric_push_level();
  }
  comm_free(&lc);

  // Partition within the node
  if (options->two_level) {
    options->two_level = 0;
    rsb(elements, nv, options, &nc, bfr);
    comm_free(&nc);
  }

  check_rsb_partition(gc, options);

  return 0;
}
