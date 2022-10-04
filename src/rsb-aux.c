#include "genmap-impl.h"
#include "metrics.h"
#include "sort.h"
#include <limits.h>
#include <time.h>

extern int rcb(struct array *elements, size_t unit_size, int ndim,
               struct comm *c, buffer *bfr);
extern int fiedler(struct array *elements, int nv, parrsb_options *options,
                   struct comm *gsc, buffer *buf, int verbose);

extern uint get_components(sint *component, struct rsb_element *elements,
                           struct comm *c, buffer *buf, uint nelt, uint nv);
extern uint get_components_v2(sint *component, uint nelt, unsigned nv,
                              struct rsb_element *elements,
                              const struct comm *ci, buffer *bfr);

// Check the bin value
static int check_bin_val(int bin, struct comm *gc) {
  if (bin < 0 || bin > 1) {
    if (gc->id == 0) {
      printf("%s:%d bin value out of range: %d\n", __FILE__, __LINE__, bin);
      fflush(stdout);
    }
    return 1;
  }
  return 0;
}

int balance_partitions(struct array *elements, int nv, struct comm *lc,
                       struct comm *gc, int bin, buffer *bfr) {
  assert(check_bin_val(bin, gc) == 0);

  struct ielem_t {
    uint index, orig;
    sint dest;
    scalar fiedler;
  };

  // Calculate expected # of elements per processor
  uint nelt = elements->n;
  slong nelgt = nelt, nglob = nelt, buf;
  comm_allreduce(lc, gs_long, gs_add, &nelgt, 1, &buf);
  comm_allreduce(gc, gs_long, gs_add, &nglob, 1, &buf);

  sint nelt_ = nglob / gc->np, nrem = nglob - nelt_ * gc->np;
  slong nelgt_exp = nelt_ * lc->np + nrem / 2 + (nrem % 2) * (1 - bin);
  slong send_cnt = nelgt - nelgt_exp > 0 ? nelgt - nelgt_exp : 0;

  // Setup gather-scatter
  uint size = nelt * nv, e, v;
  slong *ids = tcalloc(slong, size);
  struct rsb_element *elems = elements->ptr;
  for (e = 0; e < nelt; e++)
    for (v = 0; v < nv; v++)
      ids[e * nv + v] = elems[e].vertices[v];
  struct gs_data *gsh = gs_setup(ids, size, gc, 0, gs_pairwise, 0);

  sint *input = (sint *)ids;
  if (send_cnt > 0)
    for (e = 0; e < size; e++)
      input[e] = 0;
  else
    for (e = 0; e < size; e++)
      input[e] = 1;

  gs(input, gs_int, gs_add, 0, gsh, bfr);

  for (e = 0; e < nelt; e++)
    elems[e].proc = gc->id;

  sint sid = (send_cnt == 0) ? gc->id : INT_MAX, balanced = 0;
  comm_allreduce(gc, gs_int, gs_min, &sid, 1, &buf);

  struct crystal cr;

  if (send_cnt > 0) {
    struct array ielems;
    array_init(struct ielem_t, &ielems, 10);

    struct ielem_t ielem = {
        .index = 0, .orig = lc->id, .dest = -1, .fiedler = 0};
    int mul = sid == 0 ? 1 : -1;
    for (e = 0; e < nelt; e++) {
      for (v = 0; v < nv; v++) {
        if (input[e * nv + v] > 0) {
          ielem.index = e, ielem.fiedler = mul * elems[e].fiedler;
          array_cat(struct ielem_t, &ielems, &ielem, 1);
          break;
        }
      }
    }

    // Sort based on fiedler value and sets `orig` field
    parallel_sort(struct ielem_t, &ielems, fiedler, gs_double, 0, 1, lc, bfr);

    slong out[2][1], bfr[2][1], nielems = ielems.n;
    comm_scan(out, lc, gs_long, gs_add, &nielems, 1, bfr);
    slong start = out[0][0];

    sint P = gc->np - lc->np;
    sint part_size = (send_cnt + P - 1) / P;

    if (out[1][0] >= send_cnt) {
      balanced = 1;
      struct ielem_t *ptr = ielems.ptr;
      for (e = 0; start + e < send_cnt && e < ielems.n; e++)
        ptr[e].dest = sid + (start + e) / part_size;

      crystal_init(&cr, lc);
      sarray_transfer(struct ielem_t, &ielems, orig, 0, &cr);
      crystal_free(&cr);

      ptr = ielems.ptr;
      for (e = 0; e < ielems.n; e++)
        if (ptr[e].dest != -1)
          elems[ptr[e].index].proc = ptr[e].dest;
    }

    array_free(&ielems);
  }

  comm_allreduce(gc, gs_int, gs_max, &balanced, 1, &buf);
  if (balanced == 1) {
    crystal_init(&cr, gc);
    sarray_transfer(struct rsb_element, elements, proc, 0, &cr);
    crystal_free(&cr);

    // Do a load balanced sort in each partition
    parallel_sort(struct rsb_element, elements, fiedler, gs_double, 0, 1, lc,
                  bfr);
  } else {
    // Forget repair, just do a load balanced partition
    // TODO: Need to change how parallel_sort load balance
    parallel_sort(struct rsb_element, elements, fiedler, gs_double, 0, 1, gc,
                  bfr);
  }

  free(ids);
  gs_free(gsh);
}

int repair_partitions(struct array *elements, int nv, struct comm *tc,
                      struct comm *lc, int bin, struct comm *gc, buffer *bfr) {
  assert(check_bin_val(bin, gc) == 0);

  uint nelt = elements->n;

  struct rsb_element *e = elements->ptr;
  sint *cmpids = tcalloc(sint, nelt);

  int old_repair = 0;
  sint ncomp;
  if (old_repair)
    ncomp = get_components(cmpids, e, tc, bfr, nelt, nv);
  else
    ncomp = get_components_v2(cmpids, nelt, nv, e, tc, bfr);

  slong ncompg = ncomp, buf;
  comm_allreduce(lc, gs_long, gs_max, &ncompg, 1, &buf);

  sint root = (lc->id == 0) * gc->id;
  comm_allreduce(lc, gs_int, gs_max, &root, 1, &buf);

  int attempt = 0;
  int nattempts = 1;

  while (ncompg > 1 && attempt < nattempts) {
    slong *cmpcnt = tcalloc(slong, 3 * ncomp);

    uint i;
    for (i = 0; i < nelt; i++)
      cmpcnt[cmpids[i]]++;

    for (i = 0; i < ncomp; i++)
      cmpcnt[ncomp + i] = cmpcnt[i];

    comm_allreduce(tc, gs_long, gs_add, &cmpcnt[ncomp], ncomp,
                   &cmpcnt[2 * ncomp]);

    slong mincnt = LONG_MAX;
    sint minid = -1;
    for (i = 0; i < ncomp; i++) {
      if (cmpcnt[ncomp + i] < mincnt) {
        mincnt = cmpcnt[ncomp + i];
        minid = i;
      }
    }

    slong mincntg = mincnt;
    comm_allreduce(lc, gs_long, gs_min, &mincntg, 1, &buf);

    /* bin is the tie breaker */
    sint min_bin = (mincntg == mincnt) ? bin : INT_MAX;
    comm_allreduce(lc, gs_int, gs_min, &min_bin, 1, &buf);

    e = elements->ptr;
    for (i = 0; i < nelt; i++)
      e[i].proc = lc->id;

    sint low_np = (lc->np + 1) / 2;
    sint high_np = lc->np - low_np;
    sint start = (1 - bin) * low_np;
    sint P = bin * low_np + (1 - bin) * high_np;
    slong size = (mincntg + P - 1) / P;

    if (mincntg == mincnt && min_bin == bin) {
      slong in = cmpcnt[minid];
      slong out[2][1], buff[2][1];
      comm_scan(out, tc, gs_long, gs_add, &in, 1, buff);
      slong off = out[0][0];

      for (i = 0; i < nelt; i++) {
        if (cmpids[i] == minid) {
          e[i].proc = start + off / size;
          off++;
        }
      }
    }

    struct crystal cr;
    crystal_init(&cr, lc);
    sarray_transfer(struct rsb_element, elements, proc, 0, &cr);
    crystal_free(&cr);

    attempt++;

    /* Do a load balanced sort in each partition */
    parallel_sort(struct rsb_element, elements, fiedler, gs_double, 0, 1, tc,
                  bfr);

    nelt = elements->n;
    cmpids = trealloc(sint, cmpids, nelt);

    if (old_repair)
      ncomp = get_components(cmpids, elements->ptr, tc, bfr, nelt, nv);
    else
      ncomp = get_components_v2(cmpids, nelt, nv, elements->ptr, tc, bfr);

    ncompg = ncomp;
    comm_allreduce(lc, gs_long, gs_max, &ncompg, 1, &buf);

    free(cmpcnt);
  }

  free(cmpids);

  return 0;
}
static void check_rsb_partition(struct comm *gc, parrsb_options *opts) {
  int max_levels = log2ll(gc->np);
  int miter = opts->rsb_max_iter, mpass = opts->rsb_lanczos_max_restarts;

  for (int i = 0; i < max_levels; i++) {
    sint converged = 1;
    int val = (int)metric_get_value(i, RSB_FIEDLER_CALC_NITER);
    if (opts->rsb_algo == 0) {
      if (val == miter * mpass)
        converged = 0;
    } else if (opts->rsb_algo == 1) {
      if (val == mpass)
        converged = 0;
    }

    sint ibfr;
    double dbfr;
    comm_allreduce(gc, gs_int, gs_min, &converged, 1, &ibfr);
    if (converged == 0) {
      if (opts->rsb_algo == 0) {
        double final = metric_get_value(i, TOL_FNL);
        comm_allreduce(gc, gs_double, gs_min, &final, 1, &dbfr);

        double target = metric_get_value(i, TOL_TGT);
        comm_allreduce(gc, gs_double, gs_min, &target, 1, &dbfr);

        if (gc->id == 0) {
          printf("Warning: Lanczos only reached a tolerance of %lf given %lf "
                 "after %d x %d iterations in Level=%d!\n",
                 final, target, mpass, miter, i);
          fflush(stdout);
        }
      } else if (opts->rsb_algo == 1) {
        if (gc->id == 0) {
          printf("Warning: Inverse iteration didn't converge after %d "
                 "iterations in Level = %d\n",
                 mpass, i);
          fflush(stdout);
        }
      }
    }

    sint minc, maxc;
    minc = maxc = (sint)metric_get_value(i, RSB_COMPONENTS);
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

int rsb(struct array *elements, int nv, int check, parrsb_options *options,
        struct comm *gc, buffer *bfr) {
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

  unsigned ndim = (nv == 8) ? 3 : 2;
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
    fiedler(elements, nv, options, &lc, bfr, gc->id == 0);
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
    rsb(elements, nv, 0, options, &nc, bfr);
    comm_free(&nc);
  }

  if (check)
    check_rsb_partition(gc, options);

  return 0;
}
