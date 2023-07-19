#include "metrics.h"
#include "parrsb-impl.h"
#include "sort.h"

extern int fiedler(struct array *elements, int nv,
                   const parrsb_options *const options, struct comm *gsc,
                   buffer *buf, int verbose);

static void test_component_versions(struct array *elements, struct comm *lc,
                                    unsigned nv, unsigned lvl, buffer *bfr) {
  // Send elements to % P processor to create disconnected components.
  struct crystal cr;
  crystal_init(&cr, lc);

  struct rsb_element *pe = (struct rsb_element *)elements->ptr;
  for (unsigned e = 0; e < elements->n; e++)
    pe[e].proc = pe[e].globalId % lc->np;

  sarray_transfer(struct rsb_element, elements, proc, 1, &cr);

  struct comm tc0;
  int color = (lc->id < lc->np / 2);
  comm_split(lc, color, lc->id, &tc0);

  sint nc1 = get_components(NULL, elements, nv, &tc0, bfr, 0);
  sint nc2 = get_components_v2(NULL, elements, nv, &tc0, bfr, 0);
  if (nc1 != nc2) {
    if (tc0.id == 0) {
      fprintf(stderr, "Error: Level = %u SS BFS != MS BFS: %d %d\n", lvl, nc1,
              nc2);
      fflush(stderr);
    }
    exit(EXIT_FAILURE);
  }
  if (nc1 > 1) {
    if (tc0.id == 0)
      printf("Warning: Level = %u has %d disconnected components.\n", lvl, nc1);
    fflush(stdout);
  }

  comm_free(&tc0);
  sarray_transfer(struct rsb_element, elements, proc, 0, &cr);
  crystal_free(&cr);
}

static void check_rsb_partition(struct comm *gc,
                                const parrsb_options *const opts) {
  int max_levels = log2ll(gc->np);
  int miter = opts->rsb_max_iter, mpass = opts->rsb_max_passes;

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

    struct comm c;
    comm_split(gc, converged, gc->id, &c);

    slong bfr[4];
    if (converged == 0) {
      if (opts->rsb_algo == 0) {
        double init = metric_get_value(i, TOL_INIT);
        comm_allreduce(&c, gs_double, gs_min, &init, 1, (void *)bfr);

        double target = metric_get_value(i, TOL_TGT);
        comm_allreduce(&c, gs_double, gs_min, &target, 1, (void *)bfr);

        double final = metric_get_value(i, TOL_FNL);
        comm_allreduce(&c, gs_double, gs_min, &final, 1, (void *)bfr);
        if (c.id == 0) {
          printf("Warning: Lanczos reached a residual of %lf (target: %lf) "
                 "after %d x %d iterations in Level=%d!\n",
                 final, target, mpass, miter, i);
          fflush(stdout);
        }
      } else if (opts->rsb_algo == 1) {
        if (c.id == 0) {
          printf("Warning: Inverse iteration didn't converge after %d "
                 "iterations in Level = %d\n",
                 mpass, i);
          fflush(stdout);
        }
      }
    }
    comm_free(&c);

    sint minc, maxc;
    minc = maxc = (sint)metric_get_value(i, RSB_COMPONENTS_NCOMP);
    comm_allreduce(gc, gs_int, gs_min, &minc, 1, (void *)bfr);
    comm_allreduce(gc, gs_int, gs_max, &maxc, 1, (void *)bfr);

    if (maxc > 1 && gc->id == 0) {
      printf("Warning: Partition created %d/%d (min/max) disconnected "
             "components in Level=%d!\n",
             minc, maxc, i);
      fflush(stdout);
    }
  }
}

static int check_bin_val(int bin, struct comm *c) {
  if (bin < 0 || bin > 1) {
    if (c->id == 0) {
      printf("%s:%d bin value out of range: %d\n", __FILE__, __LINE__, bin);
      fflush(stdout);
    }
    return 1;
  }
  return 0;
}

static int balance_partitions(struct array *elements, int nv, struct comm *lc,
                              struct comm *gc, int bin, buffer *bfr) {
  // Return if there is only one processor.
  if (gc->np == 1)
    return 0;

  assert(check_bin_val(bin, gc) == 0);

  struct ielem_t {
    uint index, orig;
    sint dest;
    scalar fiedler;
  };

  // Calculate expected # of elements per processor.
  uint ne = elements->n;
  slong nelgt = ne, nglob = ne, wrk;
  comm_allreduce(lc, gs_long, gs_add, &nelgt, 1, &wrk);
  comm_allreduce(gc, gs_long, gs_add, &nglob, 1, &wrk);

  sint ne_ = nglob / gc->np, nrem = nglob - ne_ * gc->np;
  slong nelgt_exp = ne_ * lc->np + nrem / 2 + (nrem % 2) * (1 - bin);
  slong send_cnt = nelgt - nelgt_exp > 0 ? nelgt - nelgt_exp : 0;

  // Setup gather-scatter.
  uint size = ne * nv, e, v;
  slong *ids = tcalloc(slong, size);
  struct rsb_element *elems = (struct rsb_element *)elements->ptr;
  for (e = 0; e < ne; e++) {
    for (v = 0; v < nv; v++)
      ids[e * nv + v] = elems[e].vertices[v];
  }
  struct gs_data *gsh = gs_setup(ids, size, gc, 0, gs_pairwise, 0);

  sint *input = (sint *)ids;
  if (send_cnt > 0) {
    for (e = 0; e < size; e++)
      input[e] = 0;
  } else {
    for (e = 0; e < size; e++)
      input[e] = 1;
  }

  gs(input, gs_int, gs_add, 0, gsh, bfr);

  for (e = 0; e < ne; e++)
    elems[e].proc = gc->id;

  sint sid = (send_cnt == 0) ? gc->id : INT_MAX, balanced = 0;
  comm_allreduce(gc, gs_int, gs_min, &sid, 1, &wrk);

  struct crystal cr;

  if (send_cnt > 0) {
    struct array ielems;
    array_init(struct ielem_t, &ielems, 10);

    struct ielem_t ielem = {
        .index = 0, .orig = lc->id, .dest = -1, .fiedler = 0};
    int mul = (sid == 0) ? 1 : -1;
    for (e = 0; e < ne; e++) {
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

  comm_allreduce(gc, gs_int, gs_max, &balanced, 1, &wrk);
  if (balanced == 1) {
    crystal_init(&cr, gc);
    sarray_transfer(struct rsb_element, elements, proc, 0, &cr);
    crystal_free(&cr);

    // Do a load balanced sort in each partition
    parallel_sort(struct rsb_element, elements, fiedler, gs_double, 0, 1, lc,
                  bfr);
  } else {
    // Forget about disconnected components, just do a load balanced partition
    parallel_sort(struct rsb_element, elements, fiedler, gs_double, 0, 1, gc,
                  bfr);
  }

  free(ids), gs_free(gsh);
  return 0;
}

static int repair_partitions_v2(struct array *elems, unsigned nv,
                                struct comm *tc, struct comm *lc, unsigned bin,
                                unsigned algo, buffer *bfr) {
  assert(check_bin_val(bin, lc) == 0);

  sint nc = get_components_v2(NULL, elems, nv, tc, bfr, 0), wrk;
  comm_allreduce(lc, gs_int, gs_max, &nc, 1, &wrk);
  if (nc > 1) {
    // If nc > 1, send elements back and do RCBx, RCBy and RCBz
    struct crystal cr;
    crystal_init(&cr, lc);
    sarray_transfer(struct rsb_element, elems, proc, 0, &cr);
    crystal_free(&cr);

    // Do rcb or rib
    unsigned ndim = (nv == 8) ? 3 : 2;
    switch (algo) {
    case 0:
      parallel_sort(struct rsb_element, elems, globalId, gs_long, 0, 1, lc,
                    bfr);
      break;
    case 1:
      rcb(elems, sizeof(struct rsb_element), ndim, lc, bfr);
      break;
    case 2:
      rib(elems, sizeof(struct rsb_element), ndim, lc, bfr);
      break;
    default:
      break;
    }

    // And count number of components again. If nc > 1 still, set
    // isconnected = 1
    nc = get_components_v2(NULL, elems, nv, tc, bfr, 0);
    comm_allreduce(lc, gs_int, gs_max, &nc, 1, &wrk);
  }

  return 0;
}

static sint get_bisect_comm(struct comm *const tc, const struct comm *const lc,
                            const uint level, const uint levels,
                            const struct comm comms[3]) {
  sint pid, psize;
  if (level < levels - 1) {
    sint out[2][1], wrk[2][1], in = (comms[level + 1].id == 0);
    comm_scan(out, &comms[level], gs_int, gs_add, &in, 1, wrk);
    psize = out[1][0], pid = (comms[level + 1].id == 0) * out[0][0];
    comm_allreduce(&comms[level + 1], gs_int, gs_max, &pid, 1, wrk);
  } else {
    psize = lc->np, pid = lc->id;
  }

  const int bin = (pid >= (psize + 1) / 2);
  comm_split(lc, bin, lc->id, tc);
  return bin;
}

static uint get_level_cuts(const uint level, const uint levels,
                           const struct comm comms[3]) {
  uint n;
  if (level < levels - 1) {
    sint lvl2 = (comms[level + 1].id == 0), wrk;
    comm_allreduce(&comms[level], gs_int, gs_add, &lvl2, 1, &wrk);
    n = lvl2;
  } else {
    n = comms[level].np;
  }
  parrsb_print(&comms[0], 1, "\tLevel=%d,n=%u\n", level + 1, n);

  sint cuts = 0, pow2 = 1;
  while (pow2 < n)
    pow2 <<= 1, cuts++;

  sint wrk;
  comm_allreduce(&comms[0], gs_int, gs_max, &cuts, 1, &wrk);

  return cuts;
}

static void initialize_node_level_comm(struct comm *c,
                                       const struct comm *const gc) {
#ifdef MPI
  MPI_Comm node;
  MPI_Comm_split_type(gc->c, MPI_COMM_TYPE_SHARED, gc->id, MPI_INFO_NULL,
                      &node);
  comm_init(c, node);
  MPI_Comm_free(&node);
#else
  comm_init(1, 1);
#endif
}

void rsb(struct array *elements, int nv, const parrsb_options *const options,
         struct comm *gc, buffer *bfr) {
  // Communicator for each level of the partitioning. Right now we support a
  // maximum of two levels.inter Inter node level (level = 1) and a node level
  // (level = 2).
  struct comm comms[3];

  const int levels = options->levels;
  assert(levels >= 1 && levels <= 2);

  // Level 1 communicator is the global communicator.
  comm_dup(&comms[0], gc);
  // Node level communicator is the last level communicator.
  if (levels > 1)
    initialize_node_level_comm(&comms[levels - 1], gc);

  const sint verbose = options->verbose_level - 1;
  const uint ndim = (nv == 8) ? 3 : 2;
  for (uint level = 0; level < levels; level++) {
    // Find the maximum number of RSB cuts in current level.
    sint cuts = get_level_cuts(level, levels, comms);
    parrsb_print(gc, verbose, "Level=%u, cuts=%u\n", level + 1, cuts);

    struct comm lc;
    comm_dup(&lc, &comms[level]);
    for (uint cut = 0; cut < cuts; cut++) {
      // Run the pre-partitioner.
      parrsb_print(gc, verbose - 1, "\tPre-partition ...\n");

      metric_tic(&lc, RSB_PRE);
      switch (options->rsb_pre) {
      case 0:
        parallel_sort(struct rsb_element, elements, globalId, gs_long, 0, 1,
                      &lc, bfr);
        break;
      case 1:
        rcb(elements, sizeof(struct rsb_element), ndim, &lc, bfr);
        break;
      case 2:
        rib(elements, sizeof(struct rsb_element), ndim, &lc, bfr);
        break;
      default:
        break;
      }
      metric_toc(&lc, RSB_PRE);

      struct rsb_element *const pe = (struct rsb_element *const)elements->ptr;
      for (unsigned i = 0; i < elements->n; i++)
        pe[i].proc = lc.id;

      // Find the Fiedler vector.
      parrsb_print(gc, verbose - 1, "\tFiedler ... \n");
      metric_tic(&lc, RSB_FIEDLER);
      fiedler(elements, nv, options, &lc, bfr, gc->id == 0);
      metric_toc(&lc, RSB_FIEDLER);

      // Sort by Fiedler vector.
      parrsb_print(gc, verbose - 1, "\tSort ...\n");
      metric_tic(&lc, RSB_SORT);
      parallel_sort(struct rsb_element, elements, fiedler, gs_double, 0, 1, &lc,
                    bfr);
      metric_toc(&lc, RSB_SORT);

      // `tc` is the new communicator in newly found partitions.
      struct comm tc;
      sint bin = get_bisect_comm(&tc, &lc, level, levels, comms);

      // Find the number of disconnected components.
      parrsb_print(gc, verbose - 1, "\tComponents ...\n");
      metric_tic(&lc, RSB_COMPONENTS);
      const uint ncomp = get_components_v2(NULL, elements, nv, &tc, bfr, 0);
      metric_acc(RSB_COMPONENTS_NCOMP, ncomp);
      metric_toc(&lc, RSB_COMPONENTS);

      // Bisect and balance.
      parrsb_print(gc, verbose - 1, "\tBalance ...\n");
      metric_tic(&lc, RSB_BALANCE);
      balance_partitions(elements, nv, &tc, &lc, bin, bfr);
      metric_toc(&lc, RSB_BALANCE);

      // Split the communicator and recurse on the sub-problems.
      parrsb_print(gc, verbose - 1, "\tBisect ...\n");
      comm_free(&lc), comm_dup(&lc, &tc), comm_free(&tc);

      const uint nbrs = parrsb_get_neighbors(elements, nv, gc, &lc, bfr);
      metric_acc(RSB_NEIGHBORS, nbrs);
      metric_push_level();
    }
    comm_free(&lc);
  }
  comm_free(&comms[0]);
  if (levels > 1)
    comm_free(&comms[levels - 1]);

  check_rsb_partition(gc, options);
}
