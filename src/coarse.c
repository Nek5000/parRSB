#include "coarse.h"
#include "metrics.h"
#include <float.h>

extern int schur_setup(struct coarse *crs, struct array *eij,
                       struct crystal *cr, buffer *bfr);
extern int schur_solve(scalar *x, scalar *b, scalar tol, struct coarse *crs,
                       buffer *bfr);
extern int schur_free(struct coarse *crs);

extern void comm_split(const struct comm *old, int bin, int key,
                       struct comm *new_);

struct coarse {
  int type;
  ulong ls, lg, is, ig;
  uint ln, in, *idx;
  struct comm c;
  void *solver;
};

//------------------------------------------------------------------------------
// Number rows, local first then interface. Returns global number of local
// elements.
struct rcb_t {
  uint i, s;
  double coord[3];
  slong vtx[8];
};

static void number_local(struct array *a, uint s, uint e, int ndim, int level,
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
  int nv = (ndim == 3) ? 8 : 4;
  uint npts = size * nv;

  slong *vtx = tcalloc(slong, npts);
  for (uint i = s, k = 0; i < e; i++) {
    for (int j = 0; j < nv; j++, k++)
      vtx[k] = pa[i].vtx[j];
  }

  struct gs_data *gsh = gs_setup(vtx, npts, c, 0, gs_pairwise, 0);

  sint *dof = tcalloc(sint, npts);
  uint mid = (s + e) / 2;
  for (uint i = mid, k = (mid - s) * nv; i < e; i++) {
    for (int j = 0; j < nv; j++, k++)
      dof[k] = 1;
  }

  gs(dof, gs_int, gs_add, 0, gsh, bfr);

  for (uint i = mid, k = (mid - s) * nv; i < e; i++) {
    for (int j = 0; j < nv; j++, k++)
      dof[k] = 0;
  }

  gs(dof, gs_int, gs_add, 0, gsh, bfr);

  for (uint i = s, k = 0; i < e; i++, k++) {
    for (int j = 0; j < nv; j++) {
      if (dof[k * nv + j] > 0 && pa[i].s == INT_MAX) {
        pa[i].s = level;
        break;
      }
    }
  }

  gs_free(gsh);
  free(dof), free(vtx);

  number_local(a, s, mid, ndim, level + 1, c, bfr);
  number_local(a, mid, e, ndim, level + 1, c, bfr);
}

static void number_rows(ulong *elem, struct coarse *crs, const slong *vtx,
                        const double *coord, const uint nelt, const int nv,
                        buffer *bfr) {
  uint npts = nelt * nv;
  int nnz = (npts > 0);

  struct comm c;
  comm_split(&crs->c, nnz, crs->c.id, &c);

  uint i, j;
  if (nnz > 0) {
    int level = 1, *dof = tcalloc(int, npts);
    while (c.np > 1) {
      struct gs_data *gsh = gs_setup(vtx, npts, &c, 0, gs_pairwise, 0);

      int bin = (c.id >= (c.np + 1) / 2);
      for (i = 0; i < npts; i++)
        dof[i] = bin;

      gs(dof, gs_int, gs_add, 0, gsh, bfr);

      if (bin == 1) {
        for (i = 0; i < npts; i++)
          dof[i] = 0;
      }

      gs(dof, gs_int, gs_add, 0, gsh, bfr);

      for (i = 0; i < nelt; i++) {
        for (j = 0; j < nv; j++) {
          if (dof[i * nv + j] > 0 && !elem[i]) {
            elem[i] = level;
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

  for (i = crs->in = crs->ln = 0; i < nelt; i++) {
    if (elem[i] > 0)
      crs->in++;
    else
      crs->ln++;
  }

  slong in[2] = {crs->ln, crs->in}, out[2][2], wrk[2][2];
  comm_scan(out, &crs->c, gs_long, gs_add, in, 2, wrk);
  crs->ls = out[0][0] + 1, crs->lg = out[1][0];
  crs->is = out[0][1] + 1, crs->ig = out[1][1];

  struct array local;
  array_init(struct rcb_t, &local, crs->ln);

  int ndim = (nv == 8) ? 3 : 2;

  struct rcb_t t = {.s = INT_MAX};
  crs->idx = tcalloc(uint, nelt);
  for (uint i = 0, ln = 0, in = 0; i < nelt; i++) {
    if (elem[i] > 0)
      elem[i] = crs->lg + crs->is + in, crs->idx[crs->ln + in++] = i;
    else {
      t.i = i;
      memcpy(t.coord, &coord[i * ndim], ndim * sizeof(double));
      memcpy(t.vtx, &vtx[i * nv], nv * sizeof(slong));
      array_cat(struct rcb_t, &local, &t, 1);
      // elem[i] = crs->ls + ln, crs->idx[ln++] = i;
    }
  }

  if (local.n > 0) {
    number_local(&local, 0, local.n, ndim, 1, &c, bfr);

    sarray_sort(struct rcb_t, local.ptr, local.n, s, 0, bfr);
    struct rcb_t *pl = (struct rcb_t *)local.ptr;
    for (sint i = local.n - 1, ln = 0; i >= 0; i--)
      elem[pl[i].i] = crs->ls + ln, crs->idx[ln++] = pl[i].i;
  }

  comm_free(&c);
  array_free(&local);
}

//------------------------------------------------------------------------------
// Setup coarse grid system
//
struct coarse *coarse_setup(const unsigned int nelt, const int nv,
                            const long long *llvtx, const double *coord,
                            int type, MPI_Comm comm) {
  buffer bfr;
  buffer_init(&bfr, 1024);

  uint size = nelt * nv;
  slong *vtx = tcalloc(slong, size);
  for (uint i = 0; i < size; i++)
    vtx[i] = llvtx[i];

  struct coarse *crs = tcalloc(struct coarse, 1);
  crs->type = type;

  comm_init(&crs->c, comm);

  ulong *eid = tcalloc(ulong, nelt);
  number_rows(eid, crs, vtx, coord, nelt, nv, &bfr);

  struct crystal cr;
  crystal_init(&cr, &crs->c);

  struct array nbrs, eij;
  find_nbrs(&nbrs, eid, vtx, nelt, nv, &cr, &bfr);
  free(vtx), free(eid);
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
