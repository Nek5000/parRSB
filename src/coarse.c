#include "coarse.h"

extern int schur_setup(struct coarse *crs, struct array *eij, const ulong ng,
                       struct comm *c, struct crystal *cr, buffer *bfr);
extern int schur_solve(scalar *x, scalar *b, struct coarse *crs, struct comm *c,
                       buffer *bfr);
extern int schur_free(struct coarse *crs);

extern void comm_split(const struct comm *old, int bin, int key,
                       struct comm *new_);

//------------------------------------------------------------------------------
// Number rows, local first then interface. Returns global number of local
// elements.
static ulong number_rows(ulong *elem, ulong *ls_, ulong *lg_, uint *ln_,
                         ulong *is_, ulong *ig_, uint *in_, uint **idx_,
                         const slong *vtx, const uint nelt, const int nv,
                         const struct comm *ci, buffer *bfr) {
  uint npts = nelt * nv;
  int nnz = npts > 0;

  struct comm c;
  comm_split(ci, nnz, ci->id, &c);

  uint i, j;
  if (nnz > 0) {
    int *dof = tcalloc(int, npts), level = 1;
    while (c.np > 1) {
      struct gs_data *gsh = gs_setup(vtx, npts, &c, 0, gs_pairwise, 0);

      int bin = c.id >= (c.np + 1) / 2;
      assert(bin == 0 || bin == 1);
      for (i = 0; i < npts; i++)
        dof[i] = bin;

      gs(dof, gs_int, gs_add, 0, gsh, bfr);

      if (bin == 1)
        for (i = 0; i < npts; i++)
          dof[i] = 0;

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

  comm_free(&c);

  uint in = 0, ln = 0;
  for (i = 0; i < nelt; i++)
    if (elem[i] > 0)
      in++;
    else
      ln++;
  *in_ = in, *ln_ = ln;

  slong inp[2] = {ln, in}, out[2][2], buf[2][2];
  comm_scan(out, ci, gs_long, gs_add, inp, 2, buf);
  ulong ls = *ls_ = out[0][0] + 1, is = *is_ = out[0][1] + 1;
  ulong lg = *lg_ = out[1][0], ig = *ig_ = out[1][1];

  uint *idx = *idx_ = tcalloc(uint, nelt);
  for (i = ln = in = 0; i < nelt; i++)
    if (elem[i] > 0)
      elem[i] = lg + is++, idx[*ln_ + in++] = i;
    else
      elem[i] = ls++, idx[ln++] = i;
  assert(*ln_ == ln);      // Sanity check
  assert(*in_ == in);      // Sanity check
  assert(ln + in == nelt); // Sanity check

  return lg;
}

//------------------------------------------------------------------------------
// Setup coarse grid system
//
struct coarse *coarse_setup(uint nelt, int nv, slong const *vtx, struct comm *c,
                            buffer *bfr) {
  struct crystal cr;
  crystal_init(&cr, c);

  struct coarse *crs = tcalloc(struct coarse, 1);

  struct array nbrs, eij;
  array_init(struct nbr, &nbrs, nelt);
  array_init(struct mat_ij, &eij, nelt);

  ulong *eid = tcalloc(ulong, nelt);
  ulong ng = number_rows(eid, &crs->ls, &crs->lg, &crs->ln, &crs->is, &crs->ig,
                         &crs->in, &crs->idx, vtx, nelt, nv, c, bfr);
  find_nbrs(&nbrs, eid, vtx, nelt, nv, c, &cr, bfr);
  free(eid);

  // Convert `struct nbr` -> `struct mat_ij` and compress
  // entries which share the same (r, c) values. Set the
  // diagonal element to have zero row sum
  compress_nbrs(&eij, &nbrs, bfr);
  array_free(&nbrs);

  schur_setup(crs, &eij, ng, c, &cr, bfr);

  array_free(&eij);
  crystal_free(&cr);

  return crs;
}

int coarse_solve(scalar *x, scalar *b, struct coarse *crs, struct comm *c,
                 buffer *bfr) {
  schur_solve(x, b, crs, c, bfr);
  return 0;
}

int coarse_free(struct coarse *crs) {
  if (crs != NULL) {
    schur_free(crs);
    if (crs->idx != NULL)
      free(crs->idx), crs->idx = NULL;
    free(crs), crs = NULL;
  }
  return 0;
}
