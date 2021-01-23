#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#include <genmap-impl.h>
#include <parRSB.h>
#include <sort.h>

int genmap_rsb(genmap_handle h) {
  int verbose = h->options->debug_level > 1;
  int max_iter = 50;
  int max_pass = 50;

  genmap_comm local_c = GenmapGetLocalComm(h);
  struct comm *lc = &local_c->gsc;

  genmap_comm global_c = GenmapGetGlobalComm(h);
  struct comm *gc = &global_c->gsc;

  genmap_scan(h, local_c);
  crystal_init(&h->cr, lc);

  genmap_scan(h, GenmapGetLocalComm(h));
  uint nelt = GenmapGetNLocalElements(h);
  GenmapElements e = GenmapGetElements(h);
  GenmapInt i;
  for (i = 0; i < nelt; i++)
    e[i].globalId0 = GenmapGetLocalStartIndex(h) + i + 1;

  buffer buf;
  buffer_init(&buf, 1024);

  int nv = h->nv;
  int ndim = (nv == 8) ? 3 : 2;

  int np = gc->np;

  int level = 0, max_levels = log2(np);

  sint bfr[2];

  while (genmap_comm_size(GenmapGetLocalComm(h)) > 1) {
    local_c = GenmapGetLocalComm(h);
    lc = &local_c->gsc;
    np = lc->np;

    int global;
    if (h->options->rsb_paul == 1)
      global = 1;
    else
      global = (np == gc->np);

    /* Run RCB, RIB pre-step or just sort by global id */
    if (h->options->rsb_prepartition == 1) { // RCB
      metric_tic(lc, RCB);
      rcb(lc, h->elements, ndim, &buf);
      metric_toc(lc, RCB);
    } else if (h->options->rsb_prepartition == 2) { // RIB
      metric_tic(lc, RCB);
      rib(lc, h->elements, ndim, &buf);
      metric_toc(lc, RCB);
    } else {
      parallel_sort(struct rsb_element, h->elements, globalId0, gs_long, 0, 1,
                    lc, &buf);
    }

    /* Initialize the laplacian */
    metric_tic(lc, WEIGHTEDLAPLACIANSETUP);
    GenmapInitLaplacianWeighted(h, local_c);
    metric_toc(lc, WEIGHTEDLAPLACIANSETUP);

    /* Run fiedler */
    metric_tic(lc, FIEDLER);
    int ipass = 0, iter;
    do {
      if (h->options->rsb_algo == 0)
        iter = GenmapFiedlerLanczos(h, local_c, max_iter, global);
      else if (h->options->rsb_algo == 1)
        iter = GenmapFiedlerRQI(h, local_c, max_iter, global);
      metric_acc(NFIEDLER, iter);
      global = 0;
    } while (++ipass < max_pass && iter == max_iter);
    metric_toc(lc, FIEDLER);

    /* Dump current partition status */
    if (level < max_levels) {
      sint g_nproject, nproject;
      g_nproject = nproject = metric_get_value(level, NPROJECT);
      comm_allreduce(gc, gs_int, gs_max, &g_nproject, 1, bfr); // max

      sint g_id = (nproject == g_nproject) * gc->id;
      comm_allreduce(gc, gs_int, gs_max, &g_id, 1, bfr); // max

      sint l_id = gc->id;
      comm_allreduce(lc, gs_int, gs_max, &l_id, 1, bfr); // max

      if (g_id == l_id) {
        if (lc->id == 0)
          printf("NPROJECT=%d\n", nproject);
        // Dump the current partition
        char fname[BUFSIZ];
        sprintf(fname, "partition_level_%02d.centroids", level);
        GenmapCentroidDump(fname, h, lc->id, lc);
        sprintf(fname, "partition_level_%02d.elements", level);
        GenmapElementIdDump(fname, h, lc);

        sint discon = is_disconnected(lc, local_c->gsw, &local_c->buf, nelt, nv);
        if (discon > 0 && lc->id == 0)
          printf("\tWarning: There are disconnected components!\n");
      }
    }

    /* Sort by Fiedler vector */
    metric_tic(lc, FIEDLERSORT);
    parallel_sort(struct rsb_element, h->elements, fiedler, gs_double, 0, 1, lc,
                  &buf);
    metric_toc(lc, FIEDLERSORT);

    /* Bisect */
    metric_tic(lc, BISECT);
    int bin = 1;
    if (lc->id < (np + 1) / 2)
      bin = 0;
    // FIXME: Ugly
    GenmapSplitComm(h, &local_c, bin);
    GenmapSetLocalComm(h, local_c);
    lc = &local_c->gsc;
    genmap_scan(h, local_c);
    metric_toc(lc, BISECT);

    metric_push_level();
    level++;
  }

  /* Check if Fidler converged */
  sint converged = 1;
  for (i = 0; i < metric_get_levels(); i++) {
    int val = (int)metric_get_value(i, NFIEDLER);
    if (val >= max_pass * max_iter) {
      converged = 0;
      break;
    }
  }
  comm_allreduce(gc, gs_int, gs_min, &converged, 1, bfr); // min
  if (converged == 0 && gc->id == 0)
    printf("\tWARNING: Lanczos failed to converge while partitioning!\n");

  /* Check for disconnected components */
  GenmapInitLaplacianWeighted(h, global_c);

  sint discon = is_disconnected(gc, global_c->gsw, &global_c->buf, nelt, nv);
  if (discon > 0 && gc->id == 0)
    printf("\tWarning: There are disconnected components!\n");

  crystal_free(&h->cr);
  buffer_free(&buf);

  return 0;
}
