#include "coarse.h"
#include "genmap-impl.h"
#include "ilu.h"
#include "metrics.h"
#include "sort.h"
#include <limits.h>
#include <time.h>

extern void rcb_local(struct array *a, size_t unit_size, uint start, uint end,
                      int ndim, buffer *buf);
extern int repair_partitions(struct array *elements, int nv, struct comm *tc,
                             struct comm *lc, int bin, struct comm *gc,
                             buffer *bfr);
extern int balance_partitions(struct array *elements, int nv, struct comm *lc,
                              struct comm *gc, int bin, buffer *bfr);
extern int fiedler(struct rsb_element *elems, uint lelt, int nv, int max_iter,
                   int algo, struct comm *gsc, buffer *buf);

static void check_rsb_partition(struct comm *gc, int max_pass, int max_iter) {
  int max_levels = log2ll(gc->np);

  int i;
  for (i = 0; i < max_levels; i++) {
    sint converged = 1;
    int val = (int)metric_get_value(i, RSB_FIEDLER_NITER);
    if (val >= max_pass * max_iter)
      converged = 0;

    sint ibfr;
    comm_allreduce(gc, gs_int, gs_min, &converged, 1, &ibfr);

    if (converged == 0) {
      double dbfr;
      double final = (double)metric_get_value(i, TOL_FNL);
      comm_allreduce(gc, gs_double, gs_min, &final, 1, &dbfr);

      double target = (double)metric_get_value(i, TOL_TGT);
      comm_allreduce(gc, gs_double, gs_min, &target, 1, &dbfr);

      if (gc->id == 0) {
        printf("Warning: Partitioner only reached a tolerance of %lf given %lf "
               "after %d x %d iterations in Level=%d!\n",
               final, target, max_pass, max_iter, i);
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

static void rsb_local(struct rsb_element *elems, uint s, uint e, int nv,
                      int max_iter, int max_pass, struct comm *lc,
                      buffer *buf) {
  // lc should contain only a single rank
  if (e <= s)
    return;

  uint mid;
  uint size = e - s;
  if (size > 1) {
    fiedler(elems + s, size, nv, max_iter, 0, lc, buf);
    sarray_sort(struct rsb_element, elems + s, size, fiedler, 3, buf);
    mid = (s + e) / 2;
    rsb_local(elems, s, mid, nv, max_iter, max_pass, lc, buf);
    rsb_local(elems, mid, e, nv, max_iter, max_pass, lc, buf);
  }
}

int rsb(struct array *elements, parrsb_options *options, int nv,
        struct comm *gc, buffer *bfr) {
  int max_iter = 50;
  int max_pass = 50;

  struct comm lc;
  comm_dup(&lc, gc);

  int ndim = (nv == 8) ? 3 : 2;

  while (lc.np > 1) {
    // Run RCB, RIB pre-step or just sort by global id
    metric_tic(&lc, RSB_PRE);
    if (options->rsb_pre == 0) // Sort by global id
      parallel_sort(struct rsb_element, elements, globalId, gs_long, 0, 1, &lc,
                    bfr);
    else if (options->rsb_pre == 1) // RCB
      rcb(elements, sizeof(struct rsb_element), ndim, &lc, bfr);
    else if (options->rsb_pre == 2) // RIB
      rib(elements, sizeof(struct rsb_element), ndim, &lc, bfr);
    metric_toc(&lc, RSB_PRE);

    // Run fiedler
    metric_tic(&lc, RSB_FIEDLER);
    fiedler(elements->ptr, elements->n, nv, max_iter, options->rsb_algo, &lc,
            bfr);
    metric_toc(&lc, RSB_FIEDLER);

    // Sort by Fiedler vector
    metric_tic(&lc, RSB_FIEDLER_SORT);
    parallel_sort_2(struct rsb_element, elements, fiedler, gs_double, globalId,
                    gs_long, 0, 1, &lc, bfr);
    metric_toc(&lc, RSB_FIEDLER_SORT);

    // Bisect, repair and balance
    int bin = (lc.id >= (lc.np + 1) / 2);
    struct comm tc;
    comm_split(&lc, bin, lc.id, &tc);

    metric_tic(&lc, RSB_REPAIR_BALANCE);
    if (options->repair > 0)
      repair_partitions(elements, nv, &tc, &lc, bin, gc, bfr);
    balance_partitions(elements, nv, &tc, &lc, bin, bfr);
    metric_toc(&lc, RSB_REPAIR_BALANCE);

    comm_free(&lc);
    comm_dup(&lc, &tc);
    comm_free(&tc);

    metric_push_level();
  }

  rcb_local(elements, sizeof(struct rsb_element), 0, elements->n, 3, bfr);
  comm_free(&lc);
  check_rsb_partition(gc, max_pass, max_iter);

  return 0;
}
