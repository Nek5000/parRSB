#include "ilu.h"
#include <math.h>

struct dof {
  ulong v;
  uint d;
};

int vtx_dist(uint nelt, int nv, const slong *vtx, struct comm *c,
             struct crystal *cr, buffer *bfr) {
  uint ndofs = nelt * nv;
  int active = 0;
  if (ndofs > 0)
    active = 1;

  struct array dofs;
  array_init(struct dof, &dofs, ndofs);

  struct dof t = {0, 0};
  for (uint i = 0; i < nelt; i++) {
    for (uint j = 0; j < nv; j++) {
      t.v = vtx[i * nv + j], t.d = t.v % c->np;
      array_cat(struct dof, &dofs, &t, 1);
    }
  }

  sarray_transfer(struct dof, &dofs, d, 1, cr);
  sarray_sort_2(struct dof, dofs.ptr, dofs.n, v, 1, d, 0, bfr);

  slong count[100] = {0};
  struct dof *ptr = (struct dof *)dofs.ptr;
  if (dofs.n > 0) {
    uint np = 1;
    for (uint i = 1; i < dofs.n; i++) {
      if (ptr[i - 1].v == ptr[i].v) {
        if (ptr[i - 1].d != ptr[i].d)
          np++;
      } else
        count[np]++, np = 1;
    }
    count[np]++;
  }

  slong buf[200];
  comm_allreduce(c, gs_int, gs_add, count, 100, buf);

  if (c->id == 0) {
    for (uint i = 0; i < 100; i++)
      if (count[i] > 0)
        printf("p = %d, vtx = %lld\n", i, count[i]);
    fflush(stdout);
  }

  array_free(&dofs);

  return 0;
}

struct ev {
  ulong v;
  uint d, p;
};

struct proc {
  uint p;
};

static int elm_procs_aux(uint *p, const slong *vtx, int nv,
                         const struct array *pdst, struct array *procs,
                         buffer *bfr) {
  procs->n = 0;
  uint i, j;
  for (i = 0; i < nv; i++) {
    struct ev *ptr = (struct ev *)pdst->ptr;
    for (j = 0; j < pdst->n && ptr[j].v < vtx[i]; j++)
      ;
    for (; j < pdst->n && ptr[j].v == vtx[i]; j++)
      array_cat(struct proc, procs, &ptr[j].p, 1);
  }

  sarray_sort(struct proc, procs->ptr, procs->n, p, 0, bfr);

  struct proc *ptr = (struct proc *)procs->ptr;
  assert(procs->n > 0);
  *p = ptr[procs->n - 1].p, j = 1;
  for (i = 1; i < procs->n; i++)
    if (ptr[i - 1].p != ptr[i].p)
      j++;

  return j;
}

static int elm_procs(struct array *pdst, uint n, int nv, const slong *vtx,
                     const struct comm *c, struct crystal *cr, buffer *bfr) {
  uint ndofs = n * nv;
  int active = ndofs > 0;

  struct array dofs;
  array_init(struct dof, &dofs, ndofs);

  struct dof t = {0, 0};
  for (uint i = 0; i < n; i++) {
    for (uint j = 0; j < nv; j++) {
      t.v = vtx[i * nv + j], t.d = t.v % c->np;
      array_cat(struct dof, &dofs, &t, 1);
    }
  }

  // Make entries unique
  struct array uniq;
  array_init(struct dof, &uniq, 100);

  sarray_sort(struct dof, dofs.ptr, dofs.n, v, 1, bfr);
  struct dof *p = (struct dof *)dofs.ptr;
  uint s = 0;
  for (uint i = 1; i < dofs.n; i++) {
    if (p[s].v != p[i].v) {
      array_cat(struct dof, &uniq, &p[s], 1);
      s = i;
    }
  }
  array_cat(struct dof, &uniq, &p[s], 1);
  array_free(&dofs);

  sarray_transfer(struct dof, &uniq, d, 1, cr);
  sarray_sort_2(struct dof, uniq.ptr, uniq.n, v, 1, d, 0, bfr);

  p = (struct dof *)uniq.ptr;
  if (uniq.n > 0) {
    uint s = 0, e;
    while (s < uniq.n) {
      for (e = s + 1; e < uniq.n && p[s].v == p[e].v; e++)
        ;
      for (uint i = s; i < e; i++) {
        struct ev ev = {p[i].v, p[i].d, 0};
        for (uint j = s; j < e; j++) {
          ev.p = p[j].d;
          array_cat(struct ev, pdst, &ev, 1);
        }
      }
      s = e;
    }
  }
  array_free(&uniq);

  sarray_transfer(struct ev, pdst, d, 0, cr);
  sarray_sort_2(struct ev, pdst->ptr, pdst->n, v, 1, p, 0, bfr);
}

int elm_dist(uint n, int nv, const slong *vtx, const struct comm *c,
             struct crystal *cr, buffer *bfr) {
  struct array pdst, procs;
  array_init(struct ev, &pdst, 100);
  array_init(struct proc, &procs, 100);

  elm_procs(&pdst, n, nv, vtx, c, cr, bfr);

  slong count[100] = {0};
  uint owner;
  for (uint i = 0; i < n; i++)
    count[elm_procs_aux(&owner, &vtx[i * nv], nv, &pdst, &procs, bfr)]++;

  slong buf[100];
  comm_allreduce(c, gs_long, gs_add, count, 100, buf);

  if (c->id == 0) {
    for (uint i = 0; i < 100; i++)
      if (count[i] > 0)
        printf("p = %d, elem %lld\n", i, count[i]);
    fflush(stdout);
  }

  array_free(&procs);
  array_free(&pdst);

  return 0;
}

//=============================================================================
// ILU
//
static int local_row(const ulong *rows, const ulong row, const uint n) {
  for (uint i = 0; i < n; i++)
    if (rows[i] == row)
      return i;
  return n;
}

static int copy_row(struct array *arr, const uint i, const uint p,
                    struct par_mat *A) {
  uint *off = A->adj_off, *idx = A->adj_idx;
  ulong *rows = A->rows, *cols = A->cols;
  scalar *val = A->adj_val;

  struct mat_ij t = {rows[i], 0, 0, p, 0.0};
  for (uint j = off[i]; j < off[i + 1]; j++) {
    t.c = cols[idx[j]], t.v = val[j];
    array_cat(struct mat_ij, arr, &t, 1);
  }

  return 0;
}

static int find_owner_level(uint *owner, int *lvl, const uint n, const int nv,
                            const slong *vtx, const struct comm *c,
                            struct crystal *cr, buffer *bfr) {
  // (i, j, v)
  struct array pdst, procs;
  array_init(struct ev, &pdst, 100);
  array_init(struct proc, &procs, 100);

  elm_procs(&pdst, n, nv, vtx, c, cr, bfr);

  for (uint i = 0; i < n; i++)
    lvl[i] = elm_procs_aux(&owner[i], &vtx[i * nv], nv, &pdst, &procs, bfr);

  array_free(&procs);
  array_free(&pdst);

  return 0;
}

struct ilu {
  int nlvls, type;
  uint *lvl_off;
  struct par_mat A;
  struct crystal cr;
};

struct owner {
  ulong ri;
  uint rp, p;
};

static int ilu0_aux_rows(struct par_mat *ext, int lvl, struct ilu *ilu,
                         buffer *bfr) {
  assert(lvl > 1);
  struct par_mat *A = &ilu->A;
  assert(IS_CSR(A) && !IS_DIAG(A));

  struct crystal *cr = &ilu->cr;
  struct comm *c = &cr->comm;

  uint *off = A->adj_off, *idx = A->adj_idx, rn = A->rn;
  ulong *rows = A->rows, *cols = A->cols;
  uint *lvl_off = ilu->lvl_off;

  struct array owners, requests;
  array_init(struct owner, &owners, rn * 30);
  array_init(struct owner, &requests, rn * 30);

  struct owner t;
  for (uint i = lvl_off[lvl - 1]; i < lvl_off[lvl]; i++) {
    ulong I = rows[i];
    for (uint j = off[i]; j < off[i + 1] && cols[idx[j]] < I; j++) {
      t.ri = cols[idx[j]], t.rp = c->np, t.p = t.ri % c->np;
      array_cat(struct owner, &owners, &t, 1);
    }
  }

  for (uint i = lvl_off[0]; i < lvl_off[lvl]; i++) {
    t.ri = rows[i], t.rp = c->id, t.p = t.ri % c->np;
    array_cat(struct owner, &owners, &t, 1);
  }

  sarray_sort_2(struct owner, owners.ptr, owners.n, ri, 1, rp, 0, bfr);
  struct owner *ptr = (struct owner *)owners.ptr;
  uint i, j;
  for (i = 0; i < owners.n; i = j) {
    for (j = i + 1; j < owners.n && ptr[j].ri == ptr[i].ri; j++)
      ;
    array_cat(struct owner, &requests, &ptr[i], 1);
  }
  array_free(&owners);

  // Match row ids and set `p` to the original processor
  sarray_transfer(struct owner, &requests, p, 1, cr);

  // Set rp to the owner
  sarray_sort_2(struct owner, requests.ptr, requests.n, ri, 1, rp, 0, bfr);
  ptr = (struct owner *)requests.ptr;
  for (i = 0; i < requests.n; i = j) {
    assert(ptr[i].rp < c->np);
    for (j = i + 1; j < requests.n && ptr[j].ri == ptr[i].ri; j++) {
      assert(ptr[j].rp == c->np);
      ptr[j].rp = ptr[i].rp;
    }
  }

  // Forward requests to the owner processor
  sarray_transfer(struct owner, &requests, rp, 0, cr);

  sarray_sort_2(struct owner, requests.ptr, requests.n, ri, 1, p, 0, bfr);
  ptr = (struct owner *)requests.ptr;

  struct array sends;
  array_init(struct mat_ij, &sends, rn * 30);

  for (i = 0; i < requests.n; i = j) {
    ulong ri = ptr[i].ri;
    uint ro = local_row(rows, ri, rn);
    assert(ro < A->rn);
    for (j = i; j < requests.n && ptr[j].ri == ri; j++)
      if (ptr[j].p != c->id) // No need to send to owner
        copy_row(&sends, ro, ptr[j].p, A);
  }
  array_free(&requests);

  sarray_transfer(struct mat_ij, &sends, p, 1, cr);
  par_csr_setup(ext, &sends, 0, bfr);
  array_free(&sends);

  return 0;
}

static void ilu0_aux_update(const ulong i, const ulong k, struct par_mat *A,
                            struct par_mat *E, int verbose) {
  assert(IS_CSR(A) && !IS_DIAG(A));

  uint rn = A->rn, *off = A->adj_off, *idx = A->adj_idx, *koff, *kidx;
  scalar *val = A->adj_val, *kval;
  ulong *rows = A->rows, *cols = A->cols, *kcols;

  // Find offsets of i and k
  sint io = -1, ko = -1;
  uint j;
  for (j = 0; j < rn && (io == -1 || ko == -1); j++) {
    if (rows[j] == i)
      io = j;
    if (rows[j] == k)
      ko = j;
  }

  // Search in E if not found in A
  if (ko == -1 && E != NULL) {
    koff = E->adj_off, kidx = E->adj_idx;
    kval = E->adj_val;
    kcols = E->cols;
    ulong *krows = E->rows;
    for (j = 0; j < E->rn; j++) {
      if (krows[j] == k) {
        ko = j;
        break;
      }
    }
  } else {
    koff = off, kidx = idx;
    kval = val;
    kcols = cols;
  }

  if (io == -1 || ko == -1) {
    fprintf(stderr,
            "%s:%d ilu0_aux_update: Row not found ! i = %llu io = %d k = %llu "
            "ko = %d\n",
            __FILE__, __LINE__, i, io, k, ko);
    exit(1);
  }

  // a_ik = a_ik / a_kk
  scalar a_kk = 0;
  for (j = koff[ko]; j < koff[ko + 1]; j++) {
    if (kcols[kidx[j]] == k) {
      a_kk = kval[j];
      break;
    }
  }

  if (fabs(a_kk) < 1e-10) {
    fprintf(stderr, "%s:%d ilu0_aux_update: Diagonal is zero ! k = %llu\n",
            __FILE__, __LINE__, k);
    exit(1);
  }

  j = off[io];
  while (j < off[io + 1] && cols[idx[j]] < k)
    j++;

  if (cols[idx[j]] != k || j == off[io + 1]) {
    fprintf(stderr, "%s:%d ilu0_aux_update: a_ik is zero ! i = %llu k = %llu\n",
            __FILE__, __LINE__, i, k);
    exit(1);
  }

  scalar a_ik = val[j] / a_kk;
  if (verbose)
    printf("a_kk = %lf a_ik = %lf a_ik = %lf\n", a_kk, val[j], a_ik);
  val[j] = a_ik;

  uint kj;
  scalar a_kj;
  for (j = j + 1; j < off[io + 1]; j++) {
    for (kj = koff[ko]; kj < koff[ko + 1] && kcols[kidx[kj]] < cols[idx[j]];
         kj++)
      ;
    if (kj < koff[ko + 1] && kcols[kidx[kj]] == cols[idx[j]])
      a_kj = kval[kj];
    else
      a_kj = 0;

    if (verbose)
      printf("a_ij = %lf a_ik = %lf a_kj = %lf\n", val[j], a_ik, a_kj);
    // a_ij = a_ij - a_ik * a_kj
    val[j] -= a_ik * a_kj;
  }
}

static void ilu0_aux(int lvl, struct ilu *ilu, struct par_mat *E, int verbose,
                     buffer *bfr) {
  struct par_mat *A = &ilu->A;
  assert(IS_CSR(A) && !IS_DIAG(A));
  if (E != NULL)
    assert(IS_CSR(E) && !IS_DIAG(E));

  if (lvl < 1 || lvl > ilu->nlvls) {
    fprintf(stderr, "%s:%d ilu0_aux: level %d is not between 1 and %d\n",
            __FILE__, __LINE__, lvl, ilu->nlvls);
    exit(1);
  }

  uint *off = A->adj_off, *idx = A->adj_idx;
  scalar *val = A->adj_val;
  ulong *rows = A->rows, *cols = A->cols;

  struct crystal *cr = &ilu->cr;
  struct comm *c = &cr->comm;

  uint *lvl_off = ilu->lvl_off;
  for (uint i = lvl_off[lvl - 1] + (lvl == 1); i < lvl_off[lvl]; i++) {
    ulong I = rows[i];
    for (uint k = off[i]; k < off[i + 1] && cols[idx[k]] < I; k++) {
      ulong K = cols[idx[k]];
      if (verbose)
        printf("I = %llu K = %llu\n", I, K);
      ilu0_aux_update(I, K, A, E, verbose);
    }
  }
}

static void ilu0(struct ilu *ilu, const struct comm *c, buffer *bfr) {
  // Do ILU(0) in Level 1
  ilu0_aux(1, ilu, NULL, 0, bfr);

  for (int l = 2; l <= ilu->nlvls; l++) {
    // Ask for required rows from level i - 1
    struct par_mat E;
    ilu0_aux_rows(&E, l, ilu, bfr);
    // Perform ilu(0) at level l
    ilu0_aux(l, ilu, &E, 0, bfr);
    par_mat_free(&E);
  }

  par_mat_print(&ilu->A);
}

struct elm {
  slong vtx[8];
  uint p, lvl;
  ulong e;
};

struct ilu *ilu_setup(uint n, int nv, const slong *vtx, const struct comm *c,
                      int type, buffer *bfr) {
  assert(nv <= 8);

  struct ilu *ilu = tcalloc(struct ilu, 1);
  crystal_init(&ilu->cr, c);

  // Find owner and level for an element
  uint *owner = tcalloc(uint, n);
  int *level = tcalloc(int, n);
  find_owner_level(owner, level, n, nv, vtx, c, &ilu->cr, bfr);

  // Send the elements to owner: First do a scan to come up with the numbering
  // implied by the input and then send elements to their owner.
  slong out[2][1], buf[2][1], in = n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, buf);

  struct array elms;
  array_init(struct elm, &elms, n);

  struct elm elm;
  for (uint i = 0; i < n; i++) {
    for (int j = 0; j < nv; j++)
      elm.vtx[j] = vtx[i * nv + j];
    elm.e = out[0][0] + i + 1;
    elm.p = owner[i];
    elm.lvl = level[i];
    array_cat(struct elm, &elms, &elm, 1);
  }
  free(owner);
  free(level);

  sarray_transfer(struct elm, &elms, p, 1, &ilu->cr);
  sarray_sort_2(struct elm, elms.ptr, elms.n, lvl, 0, e, 1, bfr);

  // Setup the ILU structure: figure out global number of levels, allocate
  // data structures.
  struct elm *ptr = (struct elm *)elms.ptr;
  sint nlvls = elms.n > 0 ? ptr[elms.n - 1].lvl : 0;
  comm_allreduce(c, gs_int, gs_max, &nlvls, 1, buf);
  assert(nlvls > 0);

  ilu->type = type;
  ilu->nlvls = nlvls;
  ilu->lvl_off = (uint *)tcalloc(uint, ilu->nlvls + 1);

  ulong *eid = tcalloc(ulong, elms.n);
  uint s = 0, e = 0;
  ilu->lvl_off[0] = s;
  for (int l = 1; l <= ilu->nlvls; l++) {
    while (e < elms.n && ptr[e].lvl == l)
      e++;
    ilu->lvl_off[l] = ilu->lvl_off[l - 1] + e - s;
    s = e;
  }

  // Number rows now: All the elements in Level 0 are numbered before Level
  // 1 and so on.
  ulong ng = 0;
  for (int l = 0; l < ilu->nlvls; l++) {
    e = ilu->lvl_off[l + 1], s = ilu->lvl_off[l];
    in = e - s;
    comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
    ulong start = ng + out[0][0] + 1;
    for (; s < e; s++) {
      eid[s] = start++;
      assert(eid[s] > 0);
    }
    ng += out[1][0];
  }

  slong *vrt = tcalloc(slong, elms.n * nv);
  for (uint i = 0; i < elms.n; i++) {
    for (int j = 0; j < nv; j++) {
      vrt[i * nv + j] = ptr[i].vtx[j];
      assert(vrt[i * nv + j] > 0);
    }
  }
  array_free(&elms);

  // Find and compress neighbors in order to form the Laplacian
  struct array nbrs, eij;
  array_init(struct nbr, &nbrs, n);
  array_init(struct mat_ij, &eij, n);

  find_nbrs(&nbrs, eid, vrt, elms.n, nv, c, &ilu->cr, bfr);
  compress_nbrs(&eij, &nbrs, bfr);
  free(eid);
  free(vrt);
  array_free(&nbrs);

  // Setup the parallel CSR matrix
  par_csr_setup(&ilu->A, &eij, 0, bfr);
  array_free(&eij);

  // Setup the ILU factorization
  switch (type) {
  case 0:
    ilu0(ilu, c, bfr);
    break;
  default:
    break;
  }

  return ilu;
}

void ilu_free(struct ilu *ilu) {
  if (ilu)
    crystal_free(&ilu->cr);
  if (ilu->nlvls > 0) {
    par_mat_free(&ilu->A);
    if (ilu->lvl_off)
      free(ilu->lvl_off);
  }
  free(ilu), ilu = NULL;
}
