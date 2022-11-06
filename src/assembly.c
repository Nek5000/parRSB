#include "sort.h"

extern unsigned get_rsb_bin(uint id, uint np);

static void sfree(void *p, const char *file, unsigned line) {
  if (p)
    free(p);
}
#define tfree(p) sfree(p, __FILE__, __LINE__)

static void number_dofs_local(uint *levels, uint s, uint e, uint lvl,
                              const slong *ids, unsigned nc, struct comm *c,
                              buffer *bfr, sint *wrk) {
  if (e <= s + 1)
    return;

  uint size = e - s;
  struct gs_data *gsh = gs_setup(&ids[s], size, c, 0, gs_pairwise, 0);

  // Identify the dofs on the interface.
  uint mid = (s + e) / 2;
  for (uint i = s; i < mid; i++)
    wrk[i] = 0;
  for (uint i = mid; i < e; i++)
    wrk[i] = 1;

  gs(&wrk[s], gs_int, gs_add, 0, gsh, bfr);

  for (uint i = mid; i < e; i++)
    wrk[i] = 0;

  gs(&wrk[s], gs_int, gs_add, 0, gsh, bfr);

  for (uint i = s; i < e; i++) {
    if (wrk[i] > 0 && ids[i] > 0 && levels[i] == 0)
      levels[i] = lvl;
  }
  gs_free(gsh);

  // Recursively, go down numbering the other levels.
  number_dofs_local(levels, s, mid, lvl - 1, ids, nc, c, bfr, wrk);
  number_dofs_local(levels, mid, e, lvl - 1, ids, nc, c, bfr, wrk);
}

static void number_dofs(ulong *nid, uint n, const slong *ids, unsigned nc,
                        const struct comm *ci, buffer *bfr) {
  uint *levels = tcalloc(uint, n);
  for (uint i = 0; i < n; i++)
    levels[i] = 0;

  struct comm c;
  comm_split(ci, n > 0, ci->id, &c);

  // Let's identify the levels of the dofs as we go down the RSB partitioning
  // tree.
  if (n > 0) {
    sint *wrk = tcalloc(sint, n);
    unsigned lvl = 1e6;
    while (c.np > 1) {
      unsigned bin = get_rsb_bin(c.id, c.np);

      struct gs_data *gsh = gs_setup(ids, n, &c, 0, gs_pairwise, 0);

      for (uint i = 0; i < n; i++)
        wrk[i] = bin;
      gs(wrk, gs_int, gs_add, 0, gsh, bfr);

      if (bin) {
        for (uint i = 0; i < n; i++)
          wrk[i] = 0;
      }
      gs(wrk, gs_int, gs_add, 0, gsh, bfr);

      for (uint i = 0; i < n; i++) {
        if (wrk[i] > 0 && ids[i] > 0 && levels[i] == 0)
          levels[i] = lvl;
      }

      gs_free(gsh);

      struct comm t;
      comm_split(&c, bin, c.id, &t);
      comm_free(&c);
      comm_dup(&c, &t);
      comm_free(&t);

      lvl--;
    }
    // Now identify the levels of the local dofs. We will reduce the level to
    // 1e3 so that we can easily differentiate local and interface dofs.
    lvl = 1e3;
    number_dofs_local(levels, 0, n, lvl, ids, nc, &c, bfr, wrk);

    tfree(wrk);
  }
  comm_free(&c);

  // Number ids based on the level. We send the ids to % p to make sure the
  // numbering is consistent and continuous.
  struct dof_t {
    ulong nid, id;
    uint level, seq, p;
  };

  struct array dofs;
  array_init(struct dof_t, &dofs, n);

  struct dof_t dof;
  for (uint i = 0; i < n; i++) {
    if (ids[i]) {
      dof.id = ids[i], dof.p = dof.id % ci->np;
      dof.seq = i, dof.level = levels[i];
      array_cat(struct dof_t, &dofs, &dof, 1);
    }
  }
  tfree(levels);

  sarray_sort(struct dof_t, dofs.ptr, dofs.n, level, 0, bfr);

  sint nlvls = 0;
  if (dofs.n > 0) {
    nlvls = 1;
    struct dof_t *pd = (struct dof_t *)dofs.ptr;
    for (uint i = 1; i < dofs.n; i++)
      if (pd[i].level != pd[i - 1].level)
        nlvls++;
  }
  slong wrk[2][1];
  comm_allreduce(ci, gs_int, gs_max, &nlvls, 1, &wrk);

  struct crystal cr;
  crystal_init(&cr, ci);

  sarray_transfer(struct dof_t, &dofs, p, 1, &cr);
  sarray_sort_2(struct dof_t, dofs.ptr, dofs.n, level, 0, id, 1, bfr);

  struct dof_t *pd = (struct dof_t *)dofs.ptr;
  uint idx = 0;
  for (uint lvl = 0; lvl < nlvls; lvl++) {
    sint l = INT_MAX;
    if (idx < dofs.n)
      l = pd[idx].level;
    comm_allreduce(ci, gs_int, gs_min, &l, 1, &wrk);

    ulong id = 0;
    uint idx1 = idx, n = 0;
    for (; idx1 < dofs.n && pd[idx1].level == l; idx1++) {
      if (pd[idx1].id != id)
        id = pd[idx1].id, n++;
    }

    slong out[2][1], in = n;
    comm_scan(out, ci, gs_long, gs_add, &in, 1, wrk);
    slong s = out[0][0];

    id = 0, n = 0;
    for (; idx < idx1; idx++) {
      if (pd[idx].id != id)
        id = pd[idx].id, n++;
      pd[idx].id = s + n;
    }
  }

  sarray_transfer(struct dof_t, &dofs, p, 0, &cr);
  sarray_sort(struct dof_t, dofs.ptr, dofs.n, seq, 0, bfr);

  if (dofs.n > 0) {
    struct dof_t *pd = (struct dof_t *)dofs.ptr;
    for (uint i = 0; i < dofs.n; i++)
      nid[pd[i].seq] = pd[i].nid;
  }

  crystal_free(&cr);
  array_free(&dofs);
}

// Find compressed/unique ids (cids) and the map (u2c) from user ids (ids) to
// cids. Also, return number of compress ids (cn).
static uint unique_ids(sint *u2c, ulong *cids, uint n, const ulong *ids,
                       buffer *bfr) {
  struct id_t {
    ulong id;
    uint idx;
    sint u2c;
  };

  struct array arr;
  array_init(struct id_t, &arr, n);

  struct id_t t = {.u2c = -1};
  for (uint i = 0; i < n; i++) {
    t.id = ids[i], t.idx = i;
    array_cat(struct id_t, &arr, &t, 1);
  }

  sarray_sort(struct id_t, arr.ptr, arr.n, id, 1, bfr);

  // Ignore the ids numbered zero
  sint cn = 0;
  ulong last = 0;
  struct id_t *pa = (struct id_t *)arr.ptr;
  for (uint i = 0; i < arr.n; i++) {
    if (pa[i].id != last)
      last = cids[cn] = pa[i].id, cn++;
    pa[i].u2c = cn - 1;
  }

  sarray_sort(struct id_t, pa, n, idx, 0, bfr);
  pa = (struct id_t *)arr.ptr;
  for (uint i = 0; i < n; i++)
    u2c[i] = pa[i].u2c;

  array_free(&arr);
  return cn;
}

struct array *assembly(uint n, const ulong *ids, uint nz, const uint *Ai,
                       const uint *Aj, const double *A, unsigned null_space,
                       unsigned type, const struct comm *c) {
  comm_barrier(c);
  double t0 = comm_time();

  // May be initialize buffer with better initial size.
  buffer bfr;
  buffer_init(&bfr, 1024);

  // TODO: Do the renumbering here.
  unsigned nc = nz / n;
  ulong *nid = tcalloc(ulong, n);
  number_dofs(nid, n, ids, nc, c, &bfr);

  // TODO: Setup u2c and cid here.
  sint *u2c = tcalloc(sint, n);
  ulong *cid = tcalloc(ulong, n);
  unique_ids(u2c, cid, n, nid, &bfr);
  tfree(nid);

  struct array *mijs = tcalloc(struct array, 1);

  double t1 = comm_time() - t0, min = t1, max = t1, wrk;
  comm_allreduce(c, gs_double, gs_max, &max, 1, &wrk);
  comm_allreduce(c, gs_double, gs_min, &min, 1, &wrk);
  if (c->id == 0) {
    printf("assembly: %g %g (min max)\n", min, max);
    fflush(stdout);
  }

  tfree(u2c), tfree(cid);
  buffer_free(&bfr);

  return mijs;
}

#undef tfree
