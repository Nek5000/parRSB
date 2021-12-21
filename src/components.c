#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#include <genmap-impl.h>
#include <sort.h>

struct unmarked {
  uint index;
};

struct interface_element {
  uint index;
  uint orig;
  sint dest;
  GenmapScalar fiedler;
};

/* Check the bin value */
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

/* Find the number of disconnected components */
static sint get_components(sint *component, struct rsb_element *elements,
                           struct comm *c, buffer *buf, uint nelt, uint nv) {
  slong out[2][1], buff[2][1];
  slong nelg = nelt;
  comm_scan(out, c, gs_long, gs_add, &nelg, 1, buff);
  nelg = out[1][0];
  ulong start = out[0][0];

  if (nelg == 0)
    return 0;

  GenmapLong *p, *ids;
  GenmapMalloc(nelt * nv, &p);
  GenmapMalloc(nelt * nv, &ids);

  int null_input = 0;
  if (component == NULL) {
    GenmapMalloc(nelt, &component);
    null_input = 1;
  }

  uint e;
  for (e = 0; e < nelt; e++)
    component[e] = -1;

  struct array arr;
  array_init(struct unmarked, &arr, nelt);

  struct comm cc;

  struct unmarked u;
  uint d;
  slong nnz1, nnzg, nnzg0, nnzb;
  slong nmarked = 0;
  sint count = 0;

  do {
    /* Count unmarked elements */
    arr.n = 0;
    for (e = 0; e < nelt; e++) {
      if (component[e] == -1) {
        u.index = e;
        array_cat(struct unmarked, &arr, &u, 1);
      }
    }

    int bin = 0;
    if (arr.n > 0)
      bin = 1;
    comm_split(c, bin, c->id, &cc);

    nnz1 = nnzg = nnzg0 = 0;
    if (bin == 1) {
      /* Initialize p */
      for (e = 0; e < arr.n; e++)
        for (d = 0; d < nv; d++)
          p[e * nv + d] = 0;

      /* Mark the first non-marked element as seed */
      struct unmarked *ptr = (struct unmarked *)arr.ptr;
      slong first = start + ptr[0].index;
      slong first_ = first;
      comm_allreduce(&cc, gs_long, gs_min, &first_, 1, buff);

      if (first_ == first) {
        for (d = 0; d < nv; d++)
          p[0 * nv + d] = 1;
      }

      /* Setup gs */
      for (e = 0; e < arr.n; e++)
        for (d = 0; d < nv; d++)
          ids[e * nv + d] = elements[ptr[e].index].vertices[d];

      struct gs_data *gsh = gs_setup(ids, arr.n * nv, &cc, 0, gs_pairwise, 0);

      do {
        gs(p, gs_long, gs_add, 0, gsh, buf);

        nnz1 = 0;
        for (e = 0; e < arr.n; e++) {
          for (d = 0; d < nv; d++) {
            if (p[e * nv + d] > 0) {
              nnz1++;
              component[ptr[e].index] = count;
              break;
            }
          }

          if (d < nv) { // There was one non-zero vertex in the element
            for (d = 0; d < nv; d++)
              p[e * nv + d] = 1;
          }
        }

        nnzg0 = nnzg, nnzg = nnz1;
        comm_allreduce(&cc, gs_long, gs_add, &nnzg, 1, &nnzb);
      } while (nnzg > nnzg0);

      gs_free(gsh);
    }

    comm_free(&cc);

    comm_allreduce(c, gs_long, gs_add, &nnz1, 1, &nnzb);
    nmarked += nnz1;

    count++;
  } while (nmarked < nelg);

  GenmapFree(p);
  GenmapFree(ids);
  if (null_input == 1)
    GenmapFree(component);

  return count;
}

int balance_partitions(struct array *elements, int nv, struct comm *lc,
                       struct comm *gc, int bin, buffer *bfr) {
  assert(check_bin_val(bin, gc) == 0);

  uint nelt = elements->n;

  /* Calculate expected # of elements per processor */
  slong buf;
  slong nelgt = nelt;
  comm_allreduce(lc, gs_long, gs_add, &nelgt, 1, &buf);

  slong nglob = nelt;
  comm_allreduce(gc, gs_long, gs_add, &nglob, 1, &buf);

  slong nelt_ = nglob / gc->np;
  slong nelgt_exp = nelt_ * lc->np;
  sint nrem = nglob - nelt_ * gc->np;
  nelgt_exp += nrem / 2 + (nrem % 2) * (1 - bin);

  slong send_cnt = 0;
  if (nelgt - nelgt_exp > 0)
    send_cnt = nelgt - nelgt_exp;

  /* Setup gather-scatter */
  uint size = nelt * nv;
  slong *ids = tcalloc(slong, size);
  struct rsb_element *elems = elements->ptr;
  uint e, v;
  for (e = 0; e < nelt; e++)
    for (v = 0; v < nv; v++)
      ids[e * nv + v] = elems[e].vertices[v];

  struct gs_data *gsh = gs_setup(ids, size, gc, 0, gs_pairwise, 0);

  sint *input = tcalloc(sint, size);
  if (send_cnt > 0)
    for (e = 0; e < size; e++)
      input[e] = 0;
  else
    for (e = 0; e < size; e++)
      input[e] = 1;

  gs(input, gs_int, gs_add, 0, gsh, bfr);

  for (e = 0; e < nelt; e++)
    elems[e].proc = gc->id;

  sint start_id = (send_cnt == 0) ? gc->id : INT_MAX;
  comm_allreduce(gc, gs_int, gs_min, &start_id, 1, &buf);

  struct crystal cr;
  sint balanced = 0;

  if (send_cnt > 0) {
    int mul = -1;
    if (start_id == 0) /* we are sending to lower fiedler values */
      mul = 1;

    struct array ielems;
    array_init(struct interface_element, &ielems, 10);

    struct interface_element ielem;
    ielem.dest = -1;

    for (e = 0; e < nelt; e++) {
      for (v = 0; v < nv; v++)
        if (input[e * nv + v] > 0) {
          ielem.index = e;
          ielem.orig = lc->id;
          ielem.fiedler = mul * elems[e].fiedler;
          array_cat(struct interface_element, &ielems, &ielem, 1);
          break;
        }
    }

    parallel_sort(struct interface_element, &ielems, fiedler, gs_double, 0, 1,
                  lc, bfr);

    slong ielems_n = ielems.n;
    slong out[2][1], bfr[2][1];
    comm_scan(out, lc, gs_long, gs_add, &ielems_n, 1, bfr);
    slong start = out[0][0];

    sint P = gc->np - lc->np;
    slong part_size = (send_cnt + P - 1) / P;

    if (out[1][0] < send_cnt)
      balanced = 0;
    else {
      struct interface_element *ptr = ielems.ptr;
      for (e = 0; start + e < send_cnt && e < ielems.n; e++)
        ptr[e].dest = start_id + (start + e) / part_size;

      crystal_init(&cr, lc);
      sarray_transfer(struct interface_element, &ielems, orig, 0, &cr);
      crystal_free(&cr);

      ptr = ielems.ptr;
      for (e = 0; e < ielems.n; e++)
        if (ptr[e].dest != -1)
          elems[ptr[e].index].proc = ptr[e].dest;
    }

    array_free(&ielems);
  }

  comm_allreduce(gc, gs_int, gs_min, &balanced, 1, &buf);
  if (balanced == 1) {
    crystal_init(&cr, gc);
    sarray_transfer(struct rsb_element, elements, proc, 0, &cr);
    crystal_free(&cr);

    /* Do a load balanced sort in each partition */
    parallel_sort(struct rsb_element, elements, fiedler, gs_double, 0, 1, lc,
                  bfr);
  } else {
    /* Forget repair, just do a load balanced partition */
    parallel_sort(struct rsb_element, elements, fiedler, gs_double, 0, 1, gc,
                  bfr);
  }

  nelt = elements->n;
  sint ncomp = get_components(NULL, elements->ptr, lc, bfr, nelt, nv);
  metric_acc(COMPONENTS, ncomp);

  free(input);
  gs_free(gsh);
  free(ids);
}

int repair_partitions(struct array *elements, int nv, struct comm *tc,
                      struct comm *lc, int bin, struct comm *gc, buffer *bfr) {
  assert(check_bin_val(bin, gc) == 0);

  uint nelt = elements->n;

  struct rsb_element *e = elements->ptr;
  sint *cmpids = tcalloc(sint, nelt);
  sint ncomp = get_components(cmpids, e, tc, bfr, nelt, nv);

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
    ncomp = get_components(cmpids, elements->ptr, tc, bfr, nelt, nv);
    ncompg = ncomp;
    comm_allreduce(lc, gs_long, gs_max, &ncompg, 1, &buf);

    free(cmpcnt);
  }

  free(cmpids);

  return 0;
}