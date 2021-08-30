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

static double get_avg_nbrs(genmap_handle h, struct comm *c) {
  uint nel = h->elements->n;
  int nv = h->nv;
  uint npts = nel * nv;

  slong *ids = tcalloc(slong, npts);

  struct rsb_element *elems = h->elements->ptr;
  int e, n;
  for (e = 0; e < nel; e++) {
    for (n = 0; n < nv; n++)
      ids[e * nv + n] = elems[e].vertices[n];
  }

  struct gs_data *gsh = gs_setup(ids, npts, c, 0, gs_pairwise, 0);

  int nmsg;
  pw_data_nmsg(gsh, &nmsg);

  int b;
  comm_allreduce(c, gs_int, gs_add, &nmsg, 1, &b);

  gs_free(gsh);
  free(ids);

  return (nmsg + 1e-6) / c->np;
}

int genmap_rsb(genmap_handle h) {
  int max_iter = 50;
  int max_pass = 50;

  struct comm *lc = h->local;
  struct comm *gc = h->global;

  genmap_comm_scan(h, lc);

  uint nelt = genmap_get_nel(h);
  int nv = h->nv;
  int ndim = (nv == 8) ? 3 : 2;

  int np = gc->np;
  int level = 0;

  while ((np = lc->np) > 1) {
    int global = 1;

    /* Run RCB, RIB pre-step or just sort by global id */
    if (h->options->rsb_pre == 0) // Sort by global id
      parallel_sort(struct rsb_element, h->elements, globalId, gs_long, 0, 1,
                    lc, &h->buf);
    else if (h->options->rsb_pre == 1) // RCB
      rcb(h->elements, h->unit_size, ndim, lc, &h->buf);
    else if (h->options->rsb_pre == 2) // RIB
      rib(h->elements, h->unit_size, ndim, lc, &h->buf);

    double nbrs = get_avg_nbrs(h, lc);

    /* Initialize the laplacian */
    metric_tic(lc, LAPLACIAN_INIT);
    GenmapInitLaplacianWeighted(h, lc);
    metric_toc(lc, LAPLACIAN_INIT);

    /* Run fiedler */
    metric_tic(lc, FIEDLER);
    int ipass = 0, iter;
    do {
      if (h->options->rsb_algo == 0)
        iter = GenmapFiedlerLanczos(h, lc, max_iter, global);
      else if (h->options->rsb_algo == 1)
        iter = GenmapFiedlerRQI(h, lc, max_iter, global);
      metric_acc(FIEDLER_NITER, iter);
      global = 0;
    } while (++ipass < max_pass && iter == max_iter);
    metric_toc(lc, FIEDLER);

    /* Sort by Fiedler vector */
    parallel_sort_2(struct rsb_element, h->elements, fiedler, gs_double,
                    globalId, gs_long, 0, 1, lc, &h->buf);

    /* Bisect, repair and balance */
    int bin = 1;
    if (lc->id < (lc->np + 1) / 2)
      bin = 0;

    struct comm tc;
    genmap_comm_split(lc, bin, lc->id, &tc);
    if (h->options->repair == 1)
      repair_partitions(h, &tc, lc, bin, gc);
    balance_partitions(h, &tc, bin, lc);

    if (get_avg_nbrs(h, lc) > nbrs) {
      /* Run RCB, RIB pre-step or just sort by global id */
      if (h->options->rsb_pre == 0) // Sort by global id
        parallel_sort(struct rsb_element, h->elements, globalId, gs_long, 0, 1,
                      lc, &h->buf);
      else if (h->options->rsb_pre == 1) // RCB
        rcb(h->elements, h->unit_size, ndim, lc, &h->buf);
      else if (h->options->rsb_pre == 2) // RIB
        rib(h->elements, h->unit_size, ndim, lc, &h->buf);
    }

    comm_free(lc);
    comm_dup(lc, &tc);
    comm_free(&tc);

    genmap_comm_scan(h, lc);
    metric_push_level();
    level++;
  }

  check_rsb_partition(gc, max_pass, max_iter);

  return 0;
}
