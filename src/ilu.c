#include "ilu.h"
#include <math.h>

struct dof {
  ulong v;
  uint d;
};

int vtx_dist(uint nelt, int nv, const slong *vtx, struct crystal *cr,
             buffer *bfr) {
  uint ndofs = nelt * nv;
  int active = ndofs > 0 ? 1 : 0;
  struct array dofs;
  array_init(struct dof, &dofs, ndofs);

  struct comm *c = &cr->comm;
  struct dof t = {.v = 0, .d = 0};
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

struct ev_t {
  ulong v;
  uint d, p;
};

struct p_t {
  uint p;
};

static int elm_procs(uint *p, const slong *vtx, int nv, const struct array *v2p,
                     struct array *ps, buffer *bfr) {
  ps->n = 0;
  uint i, j;
  for (i = 0; i < nv; i++) {
    struct ev_t *ptr = (struct ev_t *)v2p->ptr;
    for (j = 0; j < v2p->n && ptr[j].v < vtx[i]; j++)
      ;
    for (; j < v2p->n && ptr[j].v == vtx[i]; j++)
      array_cat(struct p_t, ps, &ptr[j].p, 1);
  }

  assert(ps->n > 0);
  sarray_sort(struct p_t, ps->ptr, ps->n, p, 0, bfr);

  struct p_t *ptr = (struct p_t *)ps->ptr;
  *p = ptr[ps->n - 1].p, j = 1;
  for (i = 1; i < ps->n; i++)
    if (ptr[i - 1].p != ptr[i].p)
      j++;

  return j;
}

static int vtx_procs(struct array *v2p, uint n, int nv, const slong *vtx,
                     struct crystal *cr, buffer *bfr) {
  uint ndofs = n * nv;
  int active = ndofs > 0 ? 1 : 0;

  struct array dofs;
  array_init(struct dof, &dofs, ndofs);

  struct comm *c = &cr->comm;
  struct dof t = {.v = 0, .d = 0};
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
        struct ev_t ev = {.v = p[i].v, .d = p[i].d, .p = 0};
        for (uint j = s; j < e; j++) {
          ev.p = p[j].d;
          array_cat(struct ev_t, v2p, &ev, 1);
        }
      }
      s = e;
    }
  }
  array_free(&uniq);

  sarray_transfer(struct ev_t, v2p, d, 0, cr);
  sarray_sort_2(struct ev_t, v2p->ptr, v2p->n, v, 1, p, 0, bfr);
}

int elm_dist(uint n, int nv, const slong *vtx, struct crystal *cr,
             buffer *bfr) {
  struct array v2p, ps;
  array_init(struct ev_t, &v2p, 100);
  array_init(struct p_t, &ps, 100);

  vtx_procs(&v2p, n, nv, vtx, cr, bfr);

  slong count[100] = {0};
  uint owner;
  for (uint i = 0; i < n; i++)
    count[elm_procs(&owner, &vtx[i * nv], nv, &v2p, &ps, bfr)]++;

  slong buf[100];
  struct comm *c = &cr->comm;
  comm_allreduce(c, gs_long, gs_add, count, 100, buf);

  if (c->id == 0) {
    for (uint i = 0; i < 100; i++)
      if (count[i] > 0)
        printf("p = %d, elem %lld\n", i, count[i]);
    fflush(stdout);
  }

  array_free(&ps);
  array_free(&v2p);

  return 0;
}

//=============================================================================
// ILU levels
//
struct key_t {
  ulong e;
  uint p;
};

struct e2n_t {
  ulong e, n;
};

static int find_unique_nbrs(struct array *e2nm, uint n, int nv,
                            const ulong *ids, const slong *vtx,
                            struct crystal *cr, buffer *bfr) {
  struct array nbrs;
  find_nbrs(&nbrs, ids, vtx, n, nv, cr, bfr);

  array_init(struct e2n_t, e2nm, n * 10);
  if (nbrs.n > 0) {
    sarray_sort_2(struct nbr, nbrs.ptr, nbrs.n, r, 1, c, 1, bfr);
    struct nbr *pn = (struct nbr *)nbrs.ptr;

    struct e2n_t en;
    uint i, j;
    for (i = 1, j = 0; i < nbrs.n; i++) {
      if ((pn[i].r != pn[j].r) || (pn[i].c != pn[j].c)) {
        en.e = pn[j].r, en.n = pn[j].c;
        array_cat(struct e2n_t, e2nm, &en, 1);
        j = i;
      }
    }
    en.e = pn[j].r, en.n = pn[j].c;
    array_cat(struct e2n_t, e2nm, &en, 1);
    array_free(&nbrs);

    sarray_sort_2(struct e2n_t, e2nm->ptr, e2nm->n, e, 1, n, 1, bfr);
  }

  return 0;
}

static int local_dof(const ulong *rows, const ulong row, const uint n) {
  for (uint i = 0; i < n; i++)
    if (rows[i] == row)
      return i;
  return n;
}

struct request_t {
  ulong r;
  uint p, o;
};

// Fill dofs array with unique dofs found in this processr
static int update_keys(struct array *keys, struct array *nbrs, const uint ln,
                       const ulong *lids, struct crystal *cr, buffer *bfr) {
  uint i, j;
  struct array temp, rqst;
  array_init(struct request_t, &temp, nbrs->n);
  array_init(struct request_t, &rqst, nbrs->n);

  struct comm *c = &cr->comm;
  struct e2n_t *pn = (struct e2n_t *)nbrs->ptr;
  struct request_t t;
  for (i = 0; i < nbrs->n; i++) {
    t.r = pn[i].n, t.p = t.r % c->np;
    t.o = (local_dof(lids, t.r, ln) < ln);
    array_cat(struct request_t, &temp, &t, 1);
  }

  struct request_t *pt = (struct request_t *)temp.ptr;
  if (temp.n > 0) {
    sarray_sort(struct request_t, temp.ptr, temp.n, r, 1, bfr);
    for (i = 1, j = 0; i < temp.n; i++) {
      if (pt[i].r != pt[j].r) {
        array_cat(struct request_t, &rqst, &pt[j], 1);
        j = i;
      }
    }
    array_cat(struct request_t, &rqst, &pt[j], 1);
  }

  sarray_transfer(struct request_t, &rqst, p, 1, cr);
  sarray_sort_2(struct request_t, rqst.ptr, rqst.n, r, 1, o, 0, bfr);

  struct request_t *pr = (struct request_t *)rqst.ptr;
  if (rqst.n > 0) {
    for (i = 1, j = 0; i < rqst.n; i++) {
      if (pr[i].r != pr[j].r) {
        // owner for dof j, j + 1, ... i - 1 is pr[i - 1].p
        assert(pr[i - 1].o == 1);
        for (; j < i; j++)
          pr[j].o = pr[i - 1].p;
        // j = i at the end
      }
    }
    assert(pr[i - 1].o == 1);
    for (; j < i; j++)
      pr[j].o = pr[i - 1].p;
  }

  sarray_transfer(struct request_t, &rqst, o, 0, cr);
  sarray_sort_2(struct request_t, rqst.ptr, rqst.n, r, 1, p, 0, bfr);

  // All the requests are forwarded correctly. Send the data back
  // to the requesting processors. Note that the requests are unique.
  struct key_t *pk = (struct key_t *)keys->ptr;
  pr = (struct request_t *)rqst.ptr;
  temp.n = 0;
  for (i = j = 0; i < rqst.n; i++) {
    while (pk[j].e < pr[i].r)
      j++;
    // Sanity check
    assert(pk[j].e == pr[i].r);
    t.o = pr[i].p;
    for (uint k = j; k < keys->n && pk[k].e == pk[j].e; k++) {
      t.r = pk[k].e, t.p = pk[k].p;
      array_cat(struct request_t, &temp, &t, 1);
    }
  }
  array_free(&rqst);

  sarray_transfer(struct request_t, &temp, o, 0, cr);
  sarray_sort_2(struct request_t, temp.ptr, temp.n, r, 1, p, 0, bfr);

  // Update the keys array. Update here is a complete rewrite.
  struct array keyt;
  array_init(struct key_t, &keyt, temp.n);

  struct key_t s;
  pt = (struct request_t *)temp.ptr;
  for (i = 0; i < ln; i++) {
    ulong e = lids[i];
    // Find `e` in the nbrs array
    for (j = 0; j < nbrs->n && pn[j].e < e; j++)
      ;
    assert(j < nbrs->n && pn[j].e == e);
    // Now go through all the neighbors and update the keys
    for (; j < nbrs->n && pn[j].e == e; j++) {
      ulong n = pn[j].n;
      // find the key of `n` in temp
      uint k = 0;
      for (; k < temp.n && pt[k].r < n; k++)
        ;
      assert(k < temp.n && pt[k].r == n);
      for (; k < temp.n && pt[k].r == n; k++) {
        s.e = e, s.p = pt[k].p;
        array_cat(struct key_t, &keyt, &s, 1);
      }
    }
  }
  array_free(&temp);

  keys->n = 0;
  if (keyt.n > 0) {
    sarray_sort_2(struct key_t, keyt.ptr, keyt.n, e, 1, p, 0, bfr);
    pk = (struct key_t *)keyt.ptr;
    for (i = 1, j = 0; i < keyt.n; i++) {
      if ((pk[i].e != pk[j].e) || (pk[i].p != pk[j].p)) {
        array_cat(struct key_t, keys, &pk[j], 1);
        j = i;
      }
    }
    array_cat(struct key_t, keys, &pk[j], 1);
  }

  array_free(&keyt);

  return 0;
}

// This routine will update `lvl_n`, `lvl_off` and `lvl_ids` with the DOF
// belongig to current level. In the process, it will remove the DOFs and their
// connectivity from ids, and vtx arrays. `n` will be adjusted to reflect
// changes.
static int find_level_01_aux(int *lvl_n, uint *lvl_off, uint *lvl_owner,
                             ulong *lvl_ids, uint *n, ulong *ids, slong *vtx,
                             int nv, struct array *keys, struct comm *c,
                             int verbose) {
  // Find the min key size locally.
  uint i, j, k;
  sint min = INT_MAX;
  struct key_t *pk = (struct key_t *)keys->ptr;
  if (keys->n > 0) {
    for (i = 1, j = 0; i < keys->n; i++) {
      if (pk[i].e != pk[j].e) {
        // Different element, update min key size if required
        min = min > i - j ? i - j : min;
        j = i;
      }
    }
  }

  sint buf[2];
  comm_allreduce(c, gs_int, gs_min, &min, 1, buf);
  if (min == INT_MAX)
    return 0;

  int lvl = *lvl_n;
  uint off = lvl_off[lvl];
  if (keys->n > 0) {
    for (i = 1, j = 0; i < keys->n; i++) {
      if (pk[i].e != pk[j].e) {
        if (i - j == min)
          lvl_ids[off] = pk[j].e, lvl_owner[off] = pk[i - 1].p, off++;
        j = i;
      }
    }
    if (i - j == min)
      lvl_ids[off] = pk[j].e, lvl_owner[off] = pk[i - 1].p, off++;
  }

  assert(lvl < 50);
  lvl++, lvl_off[lvl] = off;
  if (verbose > 1) {
    printf("id: %d |key| = %d lvl = %d size = %u\n", c->id, min, lvl,
           lvl_off[lvl] - lvl_off[lvl - 1]);
    fflush(stdout);
  }

  // Now we have to update ids and vtx. This can be done in place.
  for (i = lvl_off[lvl - 1], j = 0, k = 0; i < lvl_off[lvl]; i++, j++) {
    for (; j < *n && ids[j] < lvl_ids[i]; j++, k++) {
      ids[k] = ids[j];
      for (int v = 0; v < nv; v++)
        vtx[k * nv + v] = vtx[j * nv + v];
    }
    assert(j < *n && ids[j] == lvl_ids[i]);
  }
  for (; j < *n; j++, k++) {
    ids[k] = ids[j];
    for (int v = 0; v < nv; v++)
      vtx[k * nv + v] = vtx[j * nv + v];
  }

  *n -= lvl_off[lvl] - lvl_off[lvl - 1], *lvl_n = lvl;

  return 0;
}

static int find_levels_01(uint *lvl_off, uint *lvl_owner, ulong *lvl_ids,
                          const uint n_, const int nv, const ulong *ids_,
                          const slong *vtx_, struct crystal *cr, int verbose,
                          buffer *bfr) {
  // Copy ids and vtx since we are going to modify them
  uint n = n_;
  ulong *ids = tcalloc(ulong, n);
  slong *vtx = tcalloc(slong, n * nv);
  for (uint i = 0, j = 0; i < n; i++) {
    ids[i] = ids_[i];
    for (int v = 0; v < nv; v++, j++)
      vtx[j] = vtx_[j];
  }

  struct comm *c = &cr->comm;

  // Initialize keys: set key of each dof to the current MPI rank
  // keys array should has unique entries and should be sorted first
  // by .e and then by .p.
  struct array keys;
  array_init(struct key_t, &keys, n);
  struct key_t e2p = {.e = 0, .p = c->id};
  for (uint i = 0; i < n; i++) {
    e2p.e = ids[i];
    array_cat(struct key_t, &keys, &e2p, 1);
  }
  sarray_sort_2(struct key_t, keys.ptr, keys.n, e, 1, p, 0, bfr);
  if (verbose > 1) {
    printf("id = %d nkeys = %d\n", c->id, keys.n);
    fflush(stdout);
  }

  slong ng = n, buf[2];
  comm_allreduce(c, gs_long, gs_add, &ng, 1, buf);

  int nlvls = 0;
  struct array nbrs;
  while (ng > 0) {
    // Find unique neighbors of a DOF. DOF is a neighbor of itself.
    find_unique_nbrs(&nbrs, n, nv, ids, vtx, cr, bfr);

    // Send and receive key to/from neighbors. We forward all the requests
    // for the key of a DOF to the processor that owns the DOF and then that
    // processor takes care of the request. To do that, we first find all the
    // unique requests.
    update_keys(&keys, &nbrs, n, ids, cr, bfr);

    // Find the min key size
    // Add all the dofs with key size equal to min key size to current level
    // Update ids and vtx by removing the dofs with min key size
    find_level_01_aux(&nlvls, lvl_off, lvl_owner, lvl_ids, &n, ids, vtx, nv,
                      &keys, c, verbose);

    ng = n;
    comm_allreduce(c, gs_long, gs_add, &ng, 1, buf);
    if (verbose > 1) {
      if (c->id == 0)
        printf("lvl = %d ng = %lld\n", nlvls, ng);
      fflush(stdout);
    }
    array_free(&nbrs);
  }

  free(ids), free(vtx);

  return nlvls;
}

//=============================================================================
// ILU(0) and ILUt
//

struct ilu {
  int nlvls, type;
  uint *lvl_off;
  struct par_mat A;
  struct crystal cr;
};

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

static int ilu0_aux_rows(struct par_mat *ext, int lvl, struct ilu *ilu,
                         buffer *bfr) {
  struct owner {
    ulong ri;
    uint rp, p;
  };

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
    uint ro = local_dof(rows, ri, rn);
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

// `vtx` array is in the order of sorted element ids
static int ilu_setup_aux(struct ilu *ilu, int nlvls, uint *lvl_off,
                         uint *lvl_owner, uint *lvl_ids, const uint n,
                         const int nv, const slong *vtx, const int verbose,
                         buffer *bfr) {
  struct elm {
    slong vtx[8];
    uint p, lvl;
    ulong e;
  };

  // Send the elements in each level to the owner
  struct array elms;
  array_init(struct elm, &elms, n);

  struct elm elm;
  for (int l = 0; l < nlvls; l++) {
    for (uint i = lvl_off[l]; i < lvl_off[l + 1]; i++) {
      elm.lvl = l + 1, elm.e = lvl_ids[i], elm.p = lvl_owner[i];
      array_cat(struct elm, &elms, &elm, 1);
    }
  }
  sarray_sort(struct elm, elms.ptr, elms.n, e, 1, bfr);

  struct elm *pe = (struct elm *)elms.ptr;
  if (elms.n > 0) {
    // Sanity check
    assert(elms.n == n);
    for (uint i = 0; i < n; i++) {
      for (int v = 0; v < nv; v++)
        pe[i].vtx[v] = vtx[i * nv + v];
    }
  }

  struct crystal *cr = &ilu->cr;
  sarray_transfer(struct elm, &elms, p, 1, cr);
  sarray_sort_2(struct elm, elms.ptr, elms.n, lvl, 0, e, 1, bfr);

  // Setup the ILU structure: allocate ILU data structures.
  ilu->nlvls = nlvls;
  ilu->lvl_off = (uint *)tcalloc(uint, ilu->nlvls + 1);

  ulong *ids = trealloc(ulong, ids, elms.n);
  uint s = 0, e = 0;
  ilu->lvl_off[0] = s;
  pe = (struct elm *)elms.ptr;
  for (int l = 1; l <= ilu->nlvls; l++) {
    while (e < elms.n && pe[e].lvl == l)
      e++;
    ilu->lvl_off[l] = ilu->lvl_off[l - 1] + e - s;
    s = e;
  }

  // Number rows now: All the elements in Level 0 are numbered before Level
  // 1 and so on.
  struct comm *c = &cr->comm;
  ulong ng = 0;
  for (int l = 0; l < ilu->nlvls; l++) {
    e = ilu->lvl_off[l + 1], s = ilu->lvl_off[l];
    slong out[2][1], buf[2][1], in = e - s;
    comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
    ulong start = ng + out[0][0] + 1;
    for (; s < e; s++)
      ids[s] = start++;
    ng += out[1][0];
  }

  if (verbose > 0) {
    for (int l = 0; l < ilu->nlvls; l++) {
      printf("id = %d, lvl = %d s = %u, e = %u\n", c->id, l, ilu->lvl_off[l],
             ilu->lvl_off[l + 1]);
      fflush(stdout);
    }
  }

  slong *vrt = tcalloc(slong, elms.n * nv);
  for (uint i = 0; i < elms.n; i++) {
    for (int j = 0; j < nv; j++)
      vrt[i * nv + j] = pe[i].vtx[j];
  }
  array_free(&elms);

  // Find and compress neighbors in order to form the Laplacian
  struct array nbrs, eij;
  find_nbrs(&nbrs, ids, vrt, elms.n, nv, cr, bfr);
  free(ids), free(vrt);
  compress_nbrs(&eij, &nbrs, bfr);
  array_free(&nbrs);

  // Setup the parallel CSR matrix
  par_csr_setup(&ilu->A, &eij, 0, bfr);
  array_free(&eij);

  return 0;
}

struct ilu *ilu_setup(const uint n, const int nv, const slong *vtx,
                      const struct comm *c, const int type, const int verbose,
                      buffer *bfr) {
  struct ilu *ilu = tcalloc(struct ilu, 1);
  crystal_init(&ilu->cr, c);

  // Establish a numbering based on input
  slong out[2][1], buf[2][1], in = n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, buf);

  ulong *eid = tcalloc(ulong, n);
  for (uint i = 0; i < n; i++)
    eid[i] = out[0][0] + i + 1;

  uint *lvl_off = tcalloc(uint, 100);
  uint *lvl_owner = tcalloc(uint, n);
  ulong *lvl_ids = tcalloc(ulong, n);
  int nlvls = find_levels_01(lvl_off, lvl_owner, lvl_ids, n, nv, eid, vtx,
                             &ilu->cr, 1, bfr);
  if (verbose > 0) {
    for (int l = 0; l < nlvls; l++) {
      for (uint i = lvl_off[l]; i < lvl_off[l + 1]; i++)
        printf("id = %d e = %llu lvl = %d owner = %u\n", c->id, lvl_ids[i],
               l + 1, lvl_owner[i]);
      fflush(stdout);
    }
  }

  ilu_setup_aux(ilu, nlvls, lvl_off, lvl_owner, lvl_ids, n, nv, vtx, 1, bfr);

  // Setup the ILU factors
  ilu->type = type;
  switch (type) {
  case 0:
    ilu0(ilu, c, bfr);
    break;
  default:
    break;
  }

  free(lvl_off), free(lvl_owner), free(lvl_ids);

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
