#include "genmap-impl.h"
#include "metrics.h"
#include "sort.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

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

// Find the number of disconnected components
static uint get_components(sint *component, struct rsb_element *elements,
                           struct comm *c, buffer *buf, uint nelt, uint nv) {
  struct unmarked {
    uint index;
  };

  slong out[2][1], wrk[2][1], in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  ulong nelg = out[1][0], start = out[0][0];

  if (nelg == 0)
    return 0;

  slong *p = tcalloc(slong, nelt * nv);
  slong *ids = tcalloc(slong, nelt * nv);

  int null_input = 0;
  if (null_input = (component == NULL))
    component = tcalloc(sint, nelt);

  uint e;
  for (e = 0; e < nelt; e++)
    component[e] = -1;

  struct array arr;
  array_init(struct unmarked, &arr, nelt);

  struct comm cc;
  struct unmarked u;
  uint d, count = 0;
  slong nnz1, nnzg, nnzg0, nnzb, nmarked = 0;
  do {
    // count unmarked elements
    arr.n = 0;
    for (e = 0; e < nelt; e++) {
      if (component[e] == -1) {
        u.index = e;
        array_cat(struct unmarked, &arr, &u, 1);
      }
    }

    int bin = (arr.n > 0);
    comm_split(c, bin, c->id, &cc);

    nnz1 = nnzg = nnzg0 = 0;
    if (bin == 1) {
      // Initialize p
      for (e = 0; e < arr.n; e++)
        for (d = 0; d < nv; d++)
          p[e * nv + d] = 0;

      // Mark the first non-marked element as seed
      struct unmarked *ptr = (struct unmarked *)arr.ptr;
      slong first = start + ptr[0].index, mfirst = first;
      comm_allreduce(&cc, gs_long, gs_min, &mfirst, 1, wrk);

      if (mfirst == first) {
        for (d = 0; d < nv; d++)
          p[0 * nv + d] = 1;
      }

      // Setup gs
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
          // There was one non-zero vertex in the element
          if (d < nv) {
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
    nmarked += nnz1, count++;
  } while (nmarked < nelg);

  free(p), free(ids);
  if (null_input == 1)
    free(component);

  return count;
}

struct cmp_t {
  uint c, p, uid;
};

static sint find_or_insert(struct array *cids, struct cmp_t *t) {
  // If there are no elements in the array, insert and exit
  if (cids->n == 0) {
    array_cat(struct cmp_t, cids, t, 1);
    return -1;
  }

  // Otherwise, we will do a binary search
  struct cmp_t *pc = (struct cmp_t *)cids->ptr;
  sint s = 0, e = cids->n - 1, mid = 0;
  while (s <= e) {
    mid = (s + e) / 2;
    if (t->c == pc[mid].c)
      return pc[mid].uid;
    else if (t->c < pc[mid].c)
      e = mid - 1;
    else // t->c > pc[mid].c
      s = mid + 1;
  }

  // Okay, not found -- insert at `mid` or `mid + 1`
  uint max = cids->max;
  if (max == cids->n) {
    max += max / 2 + 1;
    pc = (struct cmp_t *)array_reserve(struct cmp_t, cids, max);
  }

  uint n = mid;
  if (t->c > pc[mid].c)
    n = mid + 1;

  struct cmp_t t0 = *t, t1;
  for (; n < cids->n; n++) {
    t1 = pc[n];
    pc[n] = t0;
    t0 = t1;
  }
  pc[n] = t0, cids->n++;

  return -1;
}

static uint get_components_v2(sint *component, uint nelt, unsigned nv,
                              struct rsb_element *elements,
                              const struct comm *ci, buffer *bfr) {
  slong out[2][1], wrk[2][1], in = nelt;
  comm_scan(out, ci, gs_long, gs_add, &in, 1, wrk);
  ulong nelg = out[1][0], start = out[0][0];

  if (nelg == 0)
    return 0;

  uint nev = nelt * nv;
  sint *p = tcalloc(sint, nev);
  slong *ids = tcalloc(slong, nev);
  uint *inds = tcalloc(uint, nev);

  int null_input = (component == NULL);
  if (null_input)
    component = tcalloc(sint, nelt);

  for (uint e = 0; e < nelt; e++)
    component[e] = -1;

  struct comm c;
  slong nmkd = 0, nc = 0;
  do {
    // Copy unmarked elements to ids
    uint unmkd = 0;
    for (uint e = 0; e < nelt; e++) {
      if (component[e] == -1) {
        assert(unmkd >= 0 && unmkd < nelt);
        inds[unmkd] = e;
        for (uint v = 0; v < nv; v++)
          ids[unmkd * nv + v] = elements[e].vertices[v];
        unmkd++;
      }
    }
    assert(unmkd <= nelt);

    int bin = (unmkd > 0);
    comm_split(ci, bin, ci->id, &c);

    slong nnzg = 0, nnzg0, ncg = 0;
    if (bin == 1) {
      // Setup gs
      struct gs_data *gsh = gs_setup(ids, unmkd * nv, &c, 0, gs_pairwise, 0);

      // Mark the first unmarked element as seed for the component c.id
      for (uint v = 0; v < nv; v++)
        p[0 * nv + v] = c.id;

      // Initialize the rest of p
      for (uint e = 1; e < unmkd; e++)
        for (uint v = 0; v < nv; v++)
          p[e * nv + v] = -1;

      do {
        nnzg0 = nnzg, nnzg = 0;

        gs(p, gs_int, gs_max, 0, gsh, bfr);

        for (uint e = 0; e < unmkd; e++) {
          uint v = 0;
          for (; v < nv; v++) {
            if (p[e * nv + v] > -1) {
              nnzg++;
              break;
            }
          }

          // There was one non-zero vertex in the element
          if (v < nv) {
            sint c = p[e * nv + v];
            for (v = 0; v < nv; v++)
              p[e * nv + v] = c;
          }
        }

        comm_allreduce(&c, gs_long, gs_add, &nnzg, 1, wrk);
      } while (nnzg0 < nnzg);
      gs_free(gsh);

      // Find unique local components and then use them to find unique
      // global components
      struct array cids;
      array_init(struct cmp_t, &cids, 100);

      struct cmp_t t;
      for (uint e = 0; e < unmkd; e++) {
        if (p[e * nv + 0] > -1) {
          t.c = p[e * nv + 0], t.p = t.c % c.np;
          find_or_insert(&cids, &t);
        }
      }
      assert(cids.n <= c.np);

      struct crystal cr;
      crystal_init(&cr, &c);
      sarray_transfer(struct cmp_t, &cids, p, 1, &cr);

      // find unique components and number them
      sarray_sort(struct cmp_t, cids.ptr, cids.n, c, 0, bfr);
      uint cnt = 0;
      if (cids.n > 0) {
        cnt++;
        struct cmp_t *pc = (struct cmp_t *)cids.ptr;
        for (uint i = 1; i < cids.n; i++) {
          if (pc[i].c > pc[i - 1].c)
            cnt++;
        }
      }

      in = cnt;
      comm_scan(out, &c, gs_long, gs_add, &in, 1, wrk);
      ulong s = out[0][0];
      ncg = out[1][0];

      if (cids.n > 0) {
        struct cmp_t *pc = (struct cmp_t *)cids.ptr;
        pc[0].uid = s;
        for (uint i = 1; i < cids.n; i++) {
          if (pc[i].c > pc[i - 1].c)
            s++;
          pc[i].uid = s;
        }
      }

      sarray_transfer(struct cmp_t, &cids, p, 0, &cr);
      crystal_free(&cr);

      sarray_sort(struct cmp_t, cids.ptr, cids.n, c, 0, bfr);
      for (uint e = 0; e < unmkd; e++) {
        if (p[e * nv + 0] > -1) {
          t.c = p[e * nv + 0];
          sint uid = find_or_insert(&cids, &t);
          assert(uid > -1);
          component[inds[e]] = nc + uid;
        }
      }

      array_free(&cids);
    }
    comm_free(&c);

    comm_allreduce(ci, gs_long, gs_max, &nnzg, 1, &wrk);
    nmkd += nnzg;
    comm_allreduce(ci, gs_long, gs_max, &ncg, 1, &wrk);
    nc += ncg;
  } while (nmkd < nelg);

  free(p), free(ids), free(inds);
  if (null_input == 1)
    free(component);

  return nc;
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
  // sint ncomp = get_components(cmpids, e, tc, bfr, nelt, nv);
  sint ncomp = get_components_v2(cmpids, nelt, nv, e, tc, bfr);

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
    // ncomp = get_components(cmpids, elements->ptr, tc, bfr, nelt, nv);
    // ncompg =  ncomp
    ncomp = get_components_v2(cmpids, nelt, nv, elements->ptr, tc, bfr);
    ncompg = ncomp;
    comm_allreduce(lc, gs_long, gs_max, &ncompg, 1, &buf);

    free(cmpcnt);
  }

  free(cmpids);

  return 0;
}
