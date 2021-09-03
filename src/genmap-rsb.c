#include <limits.h>
#include <time.h>

#include <genmap-impl.h>
#include <sort.h>

static void check_rsb_partition(struct comm *gc, int max_pass, int max_iter) {
  int max_levels = log2ll(gc->np);

  int i;
  for (i = 0; i < max_levels; i++) {
    sint converged = 1;
    int val = (int)metric_get_value(i, FIEDLER_NITER);
    if (val >= max_pass * max_iter)
      converged = 0;

    sint ibfr;
    comm_allreduce(gc, gs_int, gs_min, &converged, 1, &ibfr);

    if (converged == 0) {
      double dbfr;
      double final = (double)metric_get_value(i, LANCZOS_TOL_FINAL);
      comm_allreduce(gc, gs_double, gs_min, &final, 1, &dbfr);

      double target = (double)metric_get_value(i, LANCZOS_TOL_TARGET);
      comm_allreduce(gc, gs_double, gs_min, &target, 1, &dbfr);

      if (gc->id == 0) {
        printf("Warning: Partitioner only reached a tolerance of %lf given %lf "
               "after %d x %d iterations in Level=%d!\n",
               final, target, max_pass, max_iter, i);
        fflush(stdout);
      }
    }

    sint ncomps, minc, maxc;
    ncomps = minc = maxc = (sint)metric_get_value(i, COMPONENTS);
    comm_allreduce(gc, gs_int, gs_min, &minc, 1, &ibfr);
    comm_allreduce(gc, gs_int, gs_max, &maxc, 1, &ibfr);

    if (maxc > 1 && gc->id == 0) {
      printf("Warning: Partition created %d/%d (min/max) disconnected "
             "components.\n",
             minc, maxc);
      fflush(stdout);
    }
  }
}

static double get_avg_nbrs(struct rsb_element *elems, uint nel, int nv,
                           struct comm *c) {
  uint npts = nel * nv;
  slong *ids = tcalloc(slong, npts);

  int e, n;
  for (e = 0; e < nel; e++)
    for (n = 0; n < nv; n++)
      ids[e * nv + n] = elems[e].vertices[n];

  struct gs_data *gsh = gs_setup(ids, npts, c, 0, gs_pairwise, 0);

  int nmsg;
  pw_data_nmsg(gsh, &nmsg);

  int b;
  comm_allreduce(c, gs_int, gs_add, &nmsg, 1, &b);

  gs_free(gsh);
  free(ids);

  return (nmsg + 1e-6) / c->np;
}

static void rsb_local(struct rsb_element *elems, uint s, uint e, int nv,
                      int max_iter, int max_pass, struct comm *lc,
                      buffer *buf) {
  /* lc should contain only a single rank */
  uint mid;
  uint size = e - s;
  if (size > 1) {
    fiedler(elems + s, size, nv, max_iter, max_pass, 0, lc, buf);
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
    /* Run RCB, RIB pre-step or just sort by global id */
    if (options->rsb_pre == 0) // Sort by global id
      parallel_sort(struct rsb_element, elements, globalId, gs_long, 0, 1, &lc,
                    bfr);
    else if (options->rsb_pre == 1) // RCB
      rcb(elements, sizeof(struct rsb_element), ndim, &lc, bfr);
    else if (options->rsb_pre == 2) // RIB
      rib(elements, sizeof(struct rsb_element), ndim, &lc, bfr);

    double nbrs = get_avg_nbrs(elements->ptr, elements->n, nv, &lc);

    /* Run fiedler */
    metric_tic(&lc, FIEDLER);
    fiedler(elements->ptr, elements->n, nv, max_iter, max_pass,
            options->rsb_algo, &lc, bfr);
    metric_toc(&lc, FIEDLER);

    /* Sort by Fiedler vector */
    parallel_sort_2(struct rsb_element, elements, fiedler, gs_double, globalId,
                    gs_long, 0, 1, &lc, bfr);

    /* Bisect, repair and balance */
    int bin = 1;
    if (lc.id < (lc.np + 1) / 2)
      bin = 0;

    struct comm tc;
    genmap_comm_split(&lc, bin, lc.id, &tc);

    if (options->repair == 1) {
      repair_partitions(elements, nv, &tc, &lc, bin, gc, bfr);
      balance_partitions(elements, nv, &tc, &lc, bin, bfr);
    }

    if (get_avg_nbrs(elements->ptr, elements->n, nv, &lc) > nbrs) {
      /* Run RCB, RIB pre-step or just sort by global id */
      if (options->rsb_pre == 0) // Sort by global id
        parallel_sort(struct rsb_element, elements, globalId, gs_long, 0, 1,
                      &lc, bfr);
      else if (options->rsb_pre == 1) // RCB
        rcb(elements, sizeof(struct rsb_element), ndim, &lc, bfr);
      else if (options->rsb_pre == 2) // RIB
        rib(elements, sizeof(struct rsb_element), ndim, &lc, bfr);
    }

    comm_free(&lc);
    comm_dup(&lc, &tc);
    comm_free(&tc);

    metric_push_level();
  }
  rsb_local(elements->ptr, 0, elements->n, nv, max_iter, max_pass, &lc, bfr);
  comm_free(&lc);

  check_rsb_partition(gc, max_pass, max_iter);

  return 0;
}
