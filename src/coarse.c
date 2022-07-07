#include "coarse-impl.h"
#include "metrics.h"
#include <float.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
//------------------------------------------------------------------------------
// Setup coarse grid system. Initial dumb API.
//
// Number rows, local first then interface. Returns global number of local
// elements.
struct rcb_t {
  uint i, s;
  double coord[3];
  slong vtx[8];
};

static void nmbr_local_rcb(struct array *a, uint s, uint e, const unsigned nc,
                           const unsigned ndim, const unsigned level,
                           struct comm *c, buffer *bfr) {
  sint size = e - s;
  if (size <= 1)
    return;

  double max[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX},
         min[3] = {DBL_MAX, DBL_MAX, DBL_MAX};

  struct rcb_t *pa = (struct rcb_t *)a->ptr;
  for (uint i = s; i < e; i++) {
    for (int j = 0; j < ndim; j++) {
      if (pa[i].coord[j] < min[j])
        min[j] = pa[i].coord[j];
      if (pa[i].coord[j] > max[j])
        max[j] = pa[i].coord[j];
    }
  }

  double len = max[0] - min[0];
  int axis = 0;
  for (int j = 1; j < ndim; j++) {
    if (max[j] - min[j] > len)
      axis = j, len = max[j] - min[j];
  }

  struct rcb_t *ps = pa + s;
  switch (axis) {
  case 0:
    sarray_sort(struct rcb_t, ps, size, coord[0], 3, bfr);
    break;
  case 1:
    sarray_sort(struct rcb_t, ps, size, coord[1], 3, bfr);
    break;
  case 2:
    sarray_sort(struct rcb_t, ps, size, coord[2], 3, bfr);
    break;
  default:
    break;
  }

  // Number the elements in the interface
  uint npts = size * nc;
  slong *vtx = tcalloc(slong, npts);
  for (uint i = s, k = 0; i < e; i++) {
    for (int j = 0; j < nc; j++, k++)
      vtx[k] = pa[i].vtx[j];
  }

  struct gs_data *gsh = gs_setup(vtx, npts, c, 0, gs_pairwise, 0);

  sint *dof = tcalloc(sint, npts);
  uint mid = (s + e) / 2;
  for (uint i = mid, k = (mid - s) * nc; i < e; i++) {
    for (int j = 0; j < nc; j++, k++)
      dof[k] = 1;
  }

  gs(dof, gs_int, gs_add, 0, gsh, bfr);

  for (uint i = mid, k = (mid - s) * nc; i < e; i++) {
    for (int j = 0; j < nc; j++, k++)
      dof[k] = 0;
  }

  gs(dof, gs_int, gs_add, 0, gsh, bfr);

  for (uint i = s, k = 0; i < e; i++, k++) {
    for (int j = 0; j < nc; j++) {
      if (dof[k * nc + j] > 0 && pa[i].s == INT_MAX) {
        pa[i].s = level;
        break;
      }
    }
  }

  gs_free(gsh);
  free(dof), free(vtx);

  nmbr_local_rcb(a, s, mid, nc, ndim, level + 1, c, bfr);
  nmbr_local_rcb(a, mid, e, nc, ndim, level + 1, c, bfr);
}

// Number the DOFs internal first, faces second and all the rest (wire basket)
// next. This keeps zeros as is and renumber the positive entries in `ids`
// array.
static void number_dual_graph_dofs(ulong *dofs, struct coarse *crs, uint n,
                                   const slong *ids, uint nelt, unsigned ndim,
                                   const scalar *coord, buffer *bfr) {
  int nnz = (n > 0);
  struct comm c;
  comm_split(&crs->c, nnz, crs->c.id, &c);

  unsigned nc = n / nelt;
  uint i, j;
  if (nnz) {
    sint *dof = tcalloc(sint, n);
    int level = 1;
    while (c.np > 1) {
      struct gs_data *gsh = gs_setup(ids, n, &c, 0, gs_pairwise, 0);

      int bin = (c.id >= (c.np + 1) / 2);
      for (i = 0; i < n; i++)
        dof[i] = bin;

      gs(dof, gs_int, gs_add, 0, gsh, bfr);

      if (bin == 1) {
        for (i = 0; i < n; i++)
          dof[i] = 0;
      }

      gs(dof, gs_int, gs_add, 0, gsh, bfr);

      for (i = 0; i < nelt; i++) {
        for (j = 0; j < nc; j++) {
          if (dof[i * nc + j] > 0 && !dofs[i]) {
            dofs[i] = level;
            break;
          }
        }
      }

      gs_free(gsh);

      struct comm t;
      comm_split(&c, bin, c.id, &t);
      comm_free(&c);
      comm_dup(&c, &t);
      comm_free(&t);

      level++;
    }
    free(dof);
  }

  for (i = crs->n[0] = crs->n[1] = 0; i < nelt; i++) {
    if (dofs[i] > 0)
      crs->n[1]++;
    else
      crs->n[0]++;
  }

  slong in[2] = {crs->n[0], crs->n[1]}, out[2][2], wrk[2][2];
  comm_scan(out, &crs->c, gs_long, gs_add, in, 2, wrk);
  crs->s[0] = out[0][0] + 1, crs->ng[0] = out[1][0];
  crs->s[1] = out[0][1] + 1, crs->ng[1] = out[1][1];

  struct array local;
  array_init(struct rcb_t, &local, crs->n[0]);

  struct rcb_t t = {.s = INT_MAX};
  crs->idx = tcalloc(uint, nelt);
  for (uint i = 0, ln = 0, in = 0; i < nelt; i++) {
    if (dofs[i] > 0)
      dofs[i] = crs->ng[0] + crs->s[1] + in, crs->idx[crs->n[0] + in++] = i;
    else {
      t.i = i;
      memcpy(t.coord, &coord[i * ndim], ndim * sizeof(scalar));
      memcpy(t.vtx, &ids[i * nc], nc * sizeof(slong));
      array_cat(struct rcb_t, &local, &t, 1);
    }
  }

  if (local.n > 0) {
    nmbr_local_rcb(&local, 0, local.n, nc, ndim, 1, &c, bfr);
    sarray_sort(struct rcb_t, local.ptr, local.n, s, 0, bfr);
    struct rcb_t *pl = (struct rcb_t *)local.ptr;
    for (sint i = local.n - 1, ln = 0; i >= 0; i--)
      dofs[pl[i].i] = crs->s[0] + ln, crs->idx[ln++] = pl[i].i;
  }

  comm_free(&c);
  array_free(&local);
}

struct coarse *coarse_setup(unsigned n, unsigned nc, const long long *llvtx,
                            const scalar *coord, unsigned null_space,
                            unsigned type, struct comm *c) {
  struct coarse *crs = tcalloc(struct coarse, 1);
  crs->type = type;

  buffer bfr;
  buffer_init(&bfr, 1024);
  comm_dup(&crs->c, c);

  uint size = n * nc;
  slong *vtx = tcalloc(slong, size);
  for (uint i = 0; i < size; i++)
    vtx[i] = llvtx[i];

  ulong *dofs = tcalloc(ulong, n);
  unsigned ndim = (nc == 8) ? 3 : 2;
  number_dual_graph_dofs(dofs, crs, size, vtx, n, ndim, coord, &bfr);

  struct crystal cr;
  crystal_init(&cr, &crs->c);

  struct array nbrs, eij;
  find_nbrs(&nbrs, dofs, vtx, n, nc, &cr, &bfr);
  free(vtx), free(dofs);
  // Convert `struct nbr` -> `struct mij` and compress
  // entries which share the same (r, c) values. Set the
  // diagonal element to have zero row sum
  compress_nbrs(&eij, &nbrs, &bfr);
  array_free(&nbrs);

  switch (type) {
  case 0:
    schur_setup(crs, &eij, &cr, &bfr);
    break;
  default:
    break;
  }

  array_free(&eij);
  crystal_free(&cr);
  buffer_free(&bfr);

  return crs;
}

int coarse_solve(scalar *x, scalar *b, scalar tol, struct coarse *crs,
                 buffer *bfr) {
  metric_init();

  switch (crs->type) {
  case 0:
    schur_solve(x, b, tol, crs, bfr);
    break;
  default:
    break;
  }

  metric_push_level();
  metric_crs_print(&crs->c, 1);

  return 0;
}

int coarse_free(struct coarse *crs) {
  if (crs != NULL) {
    switch (crs->type) {
    case 0:
      schur_free(crs);
      break;
    default:
      break;
    }
    if (crs->idx != NULL)
      free(crs->idx);
    comm_free(&crs->c);
    free(crs), crs = NULL;
  }
  return 0;
}

//------------------------------------------------------------------------------
// Better API for coarse grid system.
//
static uint unique_ids(sint *perm, ulong *uid, uint n, const ulong *ids,
                       buffer *bfr) {
  struct id_t {
    ulong id;
    uint idx;
    sint perm;
  };

  struct array arr;
  array_init(struct id_t, &arr, n), arr.n = n;

  uint i;
  struct id_t *pa = (struct id_t *)arr.ptr;
  for (i = 0; i < n; i++)
    pa[i].id = ids[i], pa[i].idx = i;

  sarray_sort(struct id_t, pa, n, id, 1, bfr);
  pa = (struct id_t *)arr.ptr;

  // Ignore the ids numbered zero
  for (i = 0; i < arr.n && pa[i].id == 0; i++)
    pa[i].perm = -1;

  uint un = 0;
  ulong last = 0;
  for (; i < arr.n; i++) {
    ulong v = pa[i].id;
    if (v != last)
      last = uid[un] = v, un++;
    pa[i].perm = un - 1;
  }

  sarray_sort(struct id_t, pa, n, idx, 0, bfr);
  pa = (struct id_t *)arr.ptr;
  for (i = 0; i < n; i++)
    perm[i] = pa[i].perm;

  array_free(&arr);
  return un;
}

// Number rows, local first then interface. Returns global number of local
// elements.
struct rsb_t {
  uint i, s;
  slong vtx[8];
};

static void nmbr_local_rsb(struct array *a, uint s, uint e, const unsigned nc,
                           const unsigned level, struct comm *c, buffer *bfr) {}

static void number_dofs(ulong *dofs, struct coarse *crs, uint n,
                        const slong *ids, buffer *bfr) {
  int nnz = (n > 0);
  struct comm c;
  comm_split(&crs->c, nnz, crs->c.id, &c);

  uint i, j;
  if (nnz) {
    sint *dof = tcalloc(sint, n);
    int level = 1;
    while (c.np > 1) {
      struct gs_data *gsh = gs_setup(ids, n, &c, 0, gs_pairwise, 0);

      int bin = (c.id >= (c.np + 1) / 2);
      for (i = 0; i < n; i++)
        dof[i] = bin;

      gs(dof, gs_int, gs_add, 0, gsh, bfr);

      if (bin == 1) {
        for (i = 0; i < n; i++)
          dof[i] = 0;
      }

      gs(dof, gs_int, gs_add, 0, gsh, bfr);

      for (i = 0; i < n; i++) {
        if (dof[i] > 0 && dofs[i] == 0) {
          dofs[i] = level;
          break;
        }
      }

      gs_free(gsh);

      struct comm t;
      comm_split(&c, bin, c.id, &t);
      comm_free(&c);
      comm_dup(&c, &t);
      comm_free(&t);

      level++;
    }
    free(dof);
  }

  for (i = crs->n[0] = crs->n[1] = 0; i < n; i++) {
    if (dofs[i] > 0)
      crs->n[1]++;
    // Else dofs[i] should be zero. We want ids[i] to be non-zero as well.
    else if (ids[i] != 0)
      crs->n[0]++;
  }

  slong in[2] = {crs->n[0], crs->n[1]}, out[2][2], wrk[2][2];
  comm_scan(out, &crs->c, gs_long, gs_add, in, 2, wrk);
  crs->s[0] = out[0][0] + 1, crs->ng[0] = out[1][0];
  crs->s[1] = out[0][1] + 1, crs->ng[1] = out[1][1];

  for (uint i = 0, ln = 0, in = 0; i < n; i++) {
    // Interface DOF
    if (dofs[i] > 0)
      dofs[i] = crs->ng[0] + crs->s[1] + in++;
    // Local DOF if ids[i] != 0
    else if (ids[i] != 0)
      dofs[i] = crs->s[0] + ln++;
  }

  comm_free(&c);
}

// n  = ncr * nelt
// nz = ncr * ncr * nelt
struct coarse *crs_setup(uint n, const ulong *id, uint nz, const uint *Ai,
                         const uint *Aj, const scalar *A, unsigned ndim,
                         scalar *coord, unsigned null_space, unsigned type,
                         const struct comm *c) {
  struct coarse *crs = tcalloc(struct coarse, 1);
  crs->null_space = null_space, crs->type = type, crs->un = n;

  // Setup the buffer and duplicate the communicator.
  buffer bfr;
  buffer_init(&bfr, 1024);
  comm_dup(&crs->c, c);

  // Let's renumber the ids just in case its the schur solver. Schwarz solver
  // doesn't need re-numbering but we are going to go ahead and do it.
  slong *tid = tcalloc(slong, n);
  for (uint i = 0; i < n; i++)
    tid[i] = id[i];

  ulong *nid = tcalloc(ulong, n);
  number_dofs(nid, crs, n, tid, &bfr);

  // Find unique ids and user vector to compressed vector mapping.
  // This may not be necessary and get rid of this if not necessary.
  ulong *uid = tcalloc(ulong, n);
  crs->u2c = tcalloc(sint, n);
  crs->cn = unique_ids(crs->u2c, uid, n, nid, &bfr);

  // Now let's setup the coarse system. Create `struct mij` entries and pass
  // them into schur setup. Which processor owns the dof? All the local dofs
  // are owned by those specific preocessors -- interface dofs are owned in
  // a load balanced manner -- we will decide it when we compress the system.
  struct array mijs;
  array_init(struct mij, &mijs, n);

  struct mij m = {.r = 0, .c = 0, .idx = c->id, .p = 0, .v = 0};
  for (uint k = 0; k < nz; k++) {
    ulong i = nid[Ai[k]], j = nid[Aj[k]];
    if (i == 0 || j == 0 || A[k] == 0)
      continue;
    m.r = i, m.c = j, m.v = A[k], m.p = m.r % c->np;
    array_cat(struct mij, &mijs, &m, 1);
  }

  // Now let's assemble the matrix by sending rows with the same id to the same
  // processor. Row id `r` is sent to `r % c->np`. Then we sort first by `r`
  // then by `c`. Then we will figure out which processor owns row `r`.
  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct mij, &mijs, p, 1, &cr);

  struct array rijs;
  array_init(struct mij, &rijs, mijs.n);

  if (mijs.n > 0) {
    sarray_sort_2(struct mij, mijs.ptr, mijs.n, r, 1, c, 1, &bfr);
    struct mij *pm = (struct mij *)mijs.ptr;
    uint i, j;
    for (j = 1, i = 0; j < mijs.n; j++) {
      if ((pm[j].r != pm[i].r) || (pm[j].c != pm[i].c))
        array_cat(struct mij, &rijs, &pm[i], 1), i = j;
      else
        pm[i].v += pm[j].v, pm[i].idx = MIN(pm[i].idx, pm[j].idx);
    }
    array_cat(struct mij, &rijs, &pm[i], 1), i = j;

    // Sets the owner for each row. If the row id is less than `ng[0]` (i.e.,
    // local dof) we assign it to the processor with minimum rank id. If not,
    // we calculate the owner based on the following. We have to load balance
    // `ng[1] - ng[0]` interface dofs across `c->np` processors. We will
    // distribute the rows such that maximum load imbalance is at most 1.
    // We will have first `(c->np - nr)` processors have `ni` rows each and
    // last `nr` processors have `ni + 1` rows each.
    ulong *ng = crs->ng;
    uint ni = (ng[1] - ng[0]) / c->np, nr = ng[1] - ni * c->np, p;
    ulong ns = (c->np - nr) * ni;
    struct mij *pr = (struct mij *)rijs.ptr;
    for (j = 1, i = 0; j < rijs.n; j++) {
      if (pr[j].r != pr[i].r) {
        if (pr[i].r < ng[0])
          p = pr[i].idx;
        else {
          if (pr[i].r - ng[0] < ns)
            p = (pr[i].r - ng[0]) / ni;
          else
            p = (pr[i].r - ng[0] - ns) / (ni + 1) + c->np - nr;
        }

        for (uint k = i; k < j; k++)
          pr[k].idx = p;
      } else
        pr[i].idx = MIN(pr[i].idx, pr[j].idx);
    }
    for (uint k = i + 1; k < j; k++)
      pr[k].idx = pr[i].idx;
  }
  array_free(&mijs);

  sarray_transfer(struct mij, &rijs, idx, 1, &cr);

  crs->an = 0;
  if (rijs.n > 0) {
    sarray_sort_2(struct mij, rijs.ptr, rijs.n, r, 1, c, 1, &bfr);
    struct mij *pr = (struct mij *)rijs.ptr;
    uint j, i, an = 0, un = crs->un;
    for (j = 1, i = an = 0; j < rijs.n; j++) {
      if (pr[j].r != pr[i].r)
        uid[un + an] = -pr[i].r, an++, i = j;
    }
    uid[un + an] = -pr[i].r, an++;
    crs->an = an;
  }
  crs->c2a = gs_setup(uid, crs->un + crs->an, c, 0, gs_pairwise, 0);
  crs->n[1] = crs->an - crs->n[0];
  crs->s[1] = -uid[crs->un + crs->n[0]];

  switch (type) {
  case 0:
    schur_setup(crs, &rijs, &cr, &bfr);
    break;
  default:
    break;
  }

  array_free(&rijs), crystal_free(&cr), buffer_free(&bfr);
  free(tid), free(nid), free(uid);

  return crs;
}

void crs_solve(scalar *x, struct coarse *crs, scalar *b, scalar tol) { return; }

void crs_stats(struct coarse *crs) { return; }

void crs_free(struct coarse *crs) {
  if (crs != NULL) {
    if (crs->u2c)
      free(crs->u2c);
    gs_free(crs->c2a);
    free(crs), crs = NULL;
  }
  return;
}

#undef MIN
