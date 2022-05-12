#include "ilu.h"
#include <math.h>

extern void comm_split(const struct comm *old, int bin, int key,
                       struct comm *new_);

//=============================================================================
// ILU levels: Currently there are two methods of finding levels
//   1. Based on final element distribution among processors (dst_lvls)
//   2. Based on RSB levels while partitioning (rsb_lvls)
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
    sarray_sort_2(struct e2n_t, e2nm->ptr, e2nm->n, e, 1, n, 1, bfr);
  }
  array_free(&nbrs);

  return 0;
}

static int local_dof(const ulong *rows, const ulong row, const uint n) {
  for (uint i = 0; i < n; i++)
    if (rows[i] == row)
      return i;
  return n;
}

// Fill dofs array with unique dofs found in this processr
static int update_keys(struct array *keys, struct array *nbrs, const uint ln,
                       const ulong *lids, struct crystal *cr, buffer *bfr) {
  struct request_t {
    ulong r;
    uint p, o;
  };

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
static int dst_lvls_aux(int *lvl_n, uint *lvl_off, uint *lvl_owner,
                        ulong *lvl_ids, uint *n, ulong *ids, slong *vtx, int nv,
                        struct array *keys, struct comm *c, int verbose) {
  // Find the min key size locally.
  uint i, j, k;
  sint min = INT_MAX;
  struct key_t *pk = (struct key_t *)keys->ptr;
  if (keys->n > 0) {
    for (i = 1, j = 0; i < keys->n; i++) {
      if (pk[i].e != pk[j].e) {
        // Different element, update min key size if required
        min = (min > i - j ? i - j : min);
        j = i;
      }
    }
    min = (min > i - j ? i - j : min);
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

static int dst_lvls(uint *lvl_off, uint *lvl_owner, ulong *lvl_ids,
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

  // Initialize keys: set key of each dof to the current MPI rank.
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
    dst_lvls_aux(&nlvls, lvl_off, lvl_owner, lvl_ids, &n, ids, vtx, nv, &keys,
                 c, verbose);

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

static int rsb_lvls(uint *lvl_off, uint *lvl_owner, ulong *lvl_ids,
                    const uint n, const int nv, const ulong *ids,
                    const slong *vtx, struct comm *ci, int verbose,
                    buffer *bfr) {
  slong ng = n, buf[2];
  comm_allreduce(ci, gs_long, gs_add, &ng, 1, buf);

  // What we are going to do is identify the elements in the interface at each
  // level. These elements constitute the level of ILU. Owner of the element is
  // the processor which at least own a single vertex (possibly duplicated) of
  // the element.

  uint size = n * nv;
  sint *lvl = tcalloc(sint, 2 * n), *in = tcalloc(sint, size);
  if (lvl == NULL || in == NULL) {
    fprintf(stderr, "Failed to allocate lvl, owner and in !\n");
    exit(1);
  }

  struct comm c, t;
  comm_dup(&c, ci);

  sint nlvls = 1;
  while (c.np > 1) {
    struct gs_data *gsh = gs_setup(vtx, size, &c, 0, gs_pairwise, 0);

    int bin = (c.id >= (c.np + 1) / 2);
    for (uint i = 0; i < size; i++)
      in[i] = bin;

    gs(in, gs_int, gs_max, 0, gsh, bfr);

    if (bin == 1) {
      for (uint i = 0; i < size; i++)
        in[i] = 0;
    }

    gs(in, gs_int, gs_max, 0, gsh, bfr);

    for (uint i = 0; i < n; i++) {
      for (int j = 0; j < nv; j++) {
        if (in[i * nv + j] > 0) {
          if (lvl[i] == 0)
            lvl[i] = nlvls;
          break;
        }
      }
    }
    nlvls++;
    gs_free(gsh);

    comm_split(&c, bin, c.id, &t), comm_free(&c);
    comm_dup(&c, &t), comm_free(&t);
  }
  comm_free(&c);

  for (uint i = 0; i < n; i++) {
    if (lvl[i] == 0)
      lvl[i] = nlvls;
  }

  comm_allreduce(ci, gs_int, gs_max, &nlvls, 1, buf);
  if (verbose > 0) {
    if (ci->id == 0)
      printf("nlvls = %d\n", nlvls, ng);
  }

  // Reverse the level numbers
  for (uint i = 0; i < n; i++)
    lvl[i] = nlvls - lvl[i];

  struct gs_data *gsh = gs_setup(vtx, size, ci, 0, gs_pairwise, 0);
  sint *owner = lvl + n;
  for (int l = 0; l < nlvls; l++) {
    for (uint i = 0; i < n; i++) {
      if (lvl[i] == l)
        for (int j = 0; j < nv; j++)
          in[i * nv + j] = ci->id;
      else
        for (int j = 0; j < nv; j++)
          in[i * nv + j] = 0;
    }

    gs(in, gs_int, gs_max, 0, gsh, bfr);

    for (uint i = 0; i < n; i++) {
      if (lvl[i] == l) {
        for (int j = 0; j < nv; j++) {
          if (owner[i] < in[i * nv + j])
            owner[i] = in[i * nv + j];
        }
      }
    }
  }

  struct linfo_t {
    uint lvl, owner;
    ulong id;
  };

  struct array linfos;
  array_init(struct linfo_t, &linfos, n);

  struct linfo_t linfo = {.lvl = 0, .owner = 0, .id = 0};
  for (uint i = 0; i < n; i++) {
    linfo.lvl = lvl[i], linfo.owner = owner[i], linfo.id = ids[i];
    array_cat(struct linfo_t, &linfos, &linfo, 1);
  }
  sarray_sort(struct linfo_t, linfos.ptr, linfos.n, lvl, 0, bfr);

  if (linfos.n > 0) {
    struct linfo_t *pl = (struct linfo_t *)linfos.ptr;
    for (uint l = 0, i = 0; l < nlvls; l++) {
      for (; i < linfos.n && pl[i].lvl == l; i++)
        lvl_ids[i] = pl[i].id, lvl_owner[i] = pl[i].owner;
      lvl_off[l + 1] = i;
    }
  }

  array_free(&linfos);
  free(lvl), free(in);

  return nlvls;
}

static int find_lvls(uint *lvl_off, uint *lvl_owner, ulong *lvl_ids,
                     const uint n, const int nv, const ulong *ids,
                     const slong *vtx, int type, struct crystal *cr,
                     int verbose, buffer *bfr) {
  int nlvls = 0;
  switch (type) {
  case 0:
    nlvls = dst_lvls(lvl_off, lvl_owner, lvl_ids, n, nv, ids, vtx, cr, verbose,
                     bfr);
    break;
  case 1:
    nlvls = rsb_lvls(lvl_off, lvl_owner, lvl_ids, n, nv, ids, vtx, &cr->comm,
                     verbose, bfr);
    break;
  default:
    break;
  }
  return nlvls;
}

//=============================================================================
// Common functions for ILU(0) and ILUt
//
struct ilu {
  int nlvls, type, iter;
  uint *lvl_off;
  scalar tol;
  struct par_mat A;
  struct crystal cr;
};

inline static int copy_row(struct array *arr, const uint i, const uint p,
                           struct par_mat *A) {
  struct mat_ij t = {.r = A->rows[i], .c = 0, .idx = 0, .p = p, .v = 0.0};
  for (uint j = A->adj_off[i]; j < A->adj_off[i + 1]; j++) {
    t.c = A->cols[A->adj_idx[j]], t.v = A->adj_val[j];
    array_cat(struct mat_ij, arr, &t, 1);
  }

  return 0;
}

static int ilu_get_rows(struct par_mat *E, int lvl, struct ilu *ilu,
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

  struct array owners, requests;
  array_init(struct owner, &owners, A->rn * 30);
  array_init(struct owner, &requests, A->rn * 30);

  uint *lvl_off = ilu->lvl_off;
  struct owner t;
  for (uint i = lvl_off[lvl - 1]; i < lvl_off[lvl]; i++) {
    ulong I = A->rows[i];
    for (uint j = A->adj_off[i];
         j < A->adj_off[i + 1] && A->cols[A->adj_idx[j]] < I; j++) {
      t.ri = A->cols[A->adj_idx[j]], t.rp = c->np, t.p = t.ri % c->np;
      array_cat(struct owner, &owners, &t, 1);
    }
  }

  for (uint i = lvl_off[0]; i < lvl_off[lvl]; i++) {
    t.ri = A->rows[i], t.rp = c->id, t.p = t.ri % c->np;
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
  array_init(struct mat_ij, &sends, A->rn * 30);

  for (i = 0; i < requests.n; i = j) {
    ulong ri = ptr[i].ri;
    uint ro = local_dof(A->rows, ri, A->rn);
    assert(ro < A->rn);
    for (j = i; j < requests.n && ptr[j].ri == ri; j++) {
      // No need to send to owner
      if (ptr[j].p != c->id) {
        copy_row(&sends, ro, ptr[j].p, A);
      }
    }
  }
  array_free(&requests);

  sarray_transfer(struct mat_ij, &sends, p, 1, cr);
  par_csr_setup(E, &sends, 0, bfr);
  array_free(&sends);

  return 0;
}

//=============================================================================
// ILU(0)
//
static void ilu0_update_row(const uint io, const uint k, struct par_mat *A,
                            struct par_mat *E, int verbose, int lvl) {
  uint *off = A->adj_off, *idx = A->adj_idx;
  uint *koff = A->adj_off, *kidx = A->adj_idx;
  ulong *cols = A->cols, *kcols = A->cols;
  scalar *val = A->adj_val, *kval = A->adj_val;

  const ulong K = cols[idx[k]];
  const ulong I = A->rows[io];

  // Find offsets of K in A
  sint ko = -1;
  uint j;
  for (j = 0; j < A->rn; j++) {
    if (A->rows[j] == K) {
      ko = j;
      break;
    }
  }

  // Search in E if K is not found in A
  if (ko == -1 && E != NULL) {
    koff = E->adj_off, kidx = E->adj_idx;
    kval = E->adj_val, kcols = E->cols;
    for (j = 0; j < E->rn; j++) {
      if (E->rows[j] == K) {
        ko = j;
        break;
      }
    }
  }

  // Oops, K is no where to be found
  if (ko == -1) {
    fprintf(stderr, "%s:%d lvl = %d, k = %llu ko = %d\n", __FILE__, __LINE__,
            lvl, k, ko);
    exit(1);
  }

  // Calculate a_ik = a_ik / a_kk
  scalar a_kk = 0;
  for (j = koff[ko]; j < koff[ko + 1]; j++) {
    if (kcols[kidx[j]] == K) {
      a_kk = kval[j];
      break;
    }
  }

  if (fabs(a_kk) < 1e-10) {
    fprintf(stderr, "%s:%d ilu0: Diagonal is zero ! k = %llu\n", __FILE__,
            __LINE__, K);
    exit(1);
  }

  // cols[idx[k]] = K and val[k] = a_ik
  scalar a_ik = val[k] / a_kk;
  if (verbose) {
    printf("a_kk = %lf a_ik = %lf a_ik/a_kk = %lf\n", a_kk, val[j], a_ik);
    fflush(stdout);
  }
  val[k] = a_ik;

  uint kj;
  scalar a_kj;
  for (j = k + 1; j < off[io + 1]; j++) {
    for (kj = koff[ko]; kj < koff[ko + 1] && kcols[kidx[kj]] < cols[idx[j]];
         kj++)
      ;
    if (kj < koff[ko + 1] && kcols[kidx[kj]] == cols[idx[j]])
      a_kj = kval[kj];
    else
      a_kj = 0;

    if (verbose) {
      printf("a_ij = %lf a_ik = %lf a_kj = %lf\n", val[j], a_ik, a_kj);
      fflush(stdout);
    }
    // a_ij = a_ij - a_ik * a_kj
    val[j] -= a_ik * a_kj;
  }
}

static void ilu0_level(int lvl, uint *lvl_off, struct par_mat *A,
                       struct par_mat *E, int verbose) {
  ulong *cols = A->cols, *rows = A->rows;
  uint *off = A->adj_off, *idx = A->adj_idx, i, k;
  for (i = lvl_off[lvl - 1] + (lvl == 1); i < lvl_off[lvl]; i++)
    for (k = off[i]; k < off[i + 1] && cols[idx[k]] < rows[i]; k++)
      ilu0_update_row(i, k, A, E, verbose, lvl);
}

static void ilu0(struct ilu *ilu, const struct comm *c, buffer *bfr) {
  ilu0_level(1, ilu->lvl_off, &ilu->A, NULL, 0);
  for (int l = 2; l <= ilu->nlvls; l++) {
    struct par_mat E;
    ilu_get_rows(&E, l, ilu, bfr);
    ilu0_level(l, ilu->lvl_off, &ilu->A, &E, 0);
    par_mat_free(&E);
  }
}

//=============================================================================
// ILUt
//
static void ilut_update_row(struct array *ri, const uint io, const uint k,
                            struct par_mat *A, struct par_mat *E, int verbose,
                            buffer *bfr) {
  uint *off = A->adj_off, *idx = A->adj_idx;
  uint *koff = A->adj_off, *kidx = A->adj_idx;
  ulong *cols = A->cols, *kcols = A->cols;
  scalar *val = A->adj_val, *kval = A->adj_val;

  const ulong K = cols[idx[k]];
  const ulong I = A->rows[io];

  // Find offsets of K in A
  sint ko = -1;
  uint j, i;
  for (j = 0; j < A->rn; j++) {
    if (A->rows[j] == K) {
      ko = j;
      break;
    }
  }

  // Search in E if K is not found in A
  if (ko == -1 && E != NULL) {
    koff = E->adj_off, kidx = E->adj_idx;
    kval = E->adj_val, kcols = E->cols;
    for (j = 0; j < E->rn; j++) {
      if (E->rows[j] == K) {
        ko = j;
        break;
      }
    }
  }

  // Oops, K is no where to be found
  if (ko == -1) {
    fprintf(stderr, "%s:%d ilut: k = %llu ko = %d\n", __FILE__, __LINE__, k,
            ko);
    exit(1);
  }

  // Calculate a_ik = a_ik / a_kk
  scalar a_kk = 0;
  for (j = koff[ko]; j < koff[ko + 1]; j++) {
    if (kcols[kidx[j]] == K) {
      a_kk = kval[j];
      break;
    }
  }

  if (fabs(a_kk) < 1e-10) {
    fprintf(stderr, "%s:%d ilut: Diagonal is zero ! K = %llu\n", __FILE__,
            __LINE__, K);
    exit(1);
  }

  // We will update ri instead of val since we are keeping
  // track of the updated row I in that array.
  struct mat_ij *pr = (struct mat_ij *)ri->ptr;
  scalar a_ik = 0.0;
  for (i = 0; i < ri->n; i++) {
    if (pr[i].c == K) {
      a_ik = pr[i].v / a_kk;
      break;
    }
  }

  if (verbose) {
    printf("a_kk = %lf a_ik = %lf a_ik/a_kk = %lf\n", a_kk, pr[i].v, a_ik);
    fflush(stdout);
  }
  pr[i].v = a_ik;

  // a_ik if now calculated, let's do the a_ij = a_ij - a_ik * a_kj
  // for j > k

  // Add row K to row array, everything after the diagonal
  struct mat_ij t = {.r = I, .c = 0, .idx = 0, .p = 0, .v = 0};
  for (j = j + 1; j < koff[ko + 1]; j++) {
    t.c = kcols[kidx[j]], t.v = -a_ik * kval[j];
    array_cat(struct mat_ij, ri, &t, 1);
  }

  struct array ro;
  array_init(struct mat_ij, &ro, ri->n / 2 + 10);

  sarray_sort(struct mat_ij, ri->ptr, ri->n, c, 1, bfr);
  if (ri->n > 0) {
    pr = (struct mat_ij *)ri->ptr;
    for (i = 1, j = 0; i < ri->n; i++) {
      if (pr[i].c == pr[j].c)
        pr[j].v += pr[i].v;
      else {
        array_cat(struct mat_ij, &ro, &pr[j], 1);
        j = i;
      }
    }
    array_cat(struct mat_ij, &ro, &pr[j], 1);
  }

  ri->n = 0;
  array_cat(struct mat_ij, ri, ro.ptr, ro.n);
  array_free(&ro);
}

static void iluti_level(struct array *mij, int lvl, uint *lvl_off,
                        struct par_mat *A, struct par_mat *E, scalar tol,
                        int verbose, buffer *bfr) {
  ulong *cols = A->cols, *rows = A->rows;
  uint *off = A->adj_off, *idx = A->adj_idx, i, k, j;
  scalar *val = A->adj_val;

  struct array ri;
  array_init(struct mat_ij, &ri, 80);

  struct mat_ij t = {.r = 0, .c = 0, .idx = 0, .p = 0, .v = 0.0};

  for (i = lvl_off[lvl - 1] + (lvl == 1); i < lvl_off[lvl]; i++) {
    // Copy row I to array `ri`
    ri.n = 0, t.r = rows[i];
    for (k = off[i]; k < off[i + 1]; k++) {
      t.c = cols[idx[k]], t.v = val[k];
      array_cat(struct mat_ij, &ri, &t, 1);
    }

    // Perform ILUt
    for (k = off[i]; k < off[i + 1] && cols[idx[k]] < rows[i]; k++)
      ilut_update_row(&ri, i, k, A, E, verbose, bfr);

    struct mat_ij *pr = (struct mat_ij *)ri.ptr;
    for (k = 0; k < ri.n; k++) {
      if (fabs(pr[k].v) >= tol || pr[k].r == pr[k].c)
        array_cat(struct mat_ij, mij, &pr[k], 1);
    }

    // Update row I
    for (j = 0, k = off[i]; k < off[i + 1] && j < ri.n; j++) {
      if (pr[j].c == cols[idx[k]])
        val[k] = pr[j].v, k++;
    }
  }

  array_free(&ri);
}

static void iluti(struct ilu *ilu, const struct comm *c, buffer *bfr) {
  struct array mij;
  array_init(struct mat_ij, &mij, ilu->A.rn * 30);

  for (int i = 0; i < ilu->iter; i++) {
    struct par_mat *A = &ilu->A;
    ulong *cols = A->cols, *rows = A->rows;
    uint *off = A->adj_off, *idx = A->adj_idx, *lvl_off = ilu->lvl_off;
    scalar *val = A->adj_val;

    // if lvl == 1, add the first row to mij
    mij.n = 0;
    struct mat_ij t = {.r = 0, .c = 0, .idx = 0, .p = 0, .v = 0.0};
    if (lvl_off[0] < lvl_off[1]) {
      uint j = lvl_off[0];
      t.r = rows[j];
      for (uint k = off[j]; k < off[j + 1]; k++) {
        if (fabs(val[k]) >= ilu->tol || t.r == cols[idx[k]]) {
          t.c = cols[idx[k]], t.v = val[k];
          array_cat(struct mat_ij, &mij, &t, 1);
        }
      }
    }

    iluti_level(&mij, 1, lvl_off, A, NULL, ilu->tol, 0, bfr);
    for (int l = 2; l <= ilu->nlvls; l++) {
      struct par_mat E;
      ilu_get_rows(&E, l, ilu, bfr);
      iluti_level(&mij, l, lvl_off, A, &E, ilu->tol, 0, bfr);
      par_mat_free(&E);
    }

    par_mat_free(A);
    par_csr_setup(A, &mij, 0, bfr);
  }

  array_free(&mij);
}

//=============================================================================
// ILU API related functions
//
// `vtx` array is in the order of sorted element ids
static int ilu_setup_aux(struct ilu *ilu, int nlvls, uint *lvl_off,
                         uint *lvl_owner, ulong *lvl_ids, const uint n,
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
  ulong *ids = trealloc(ulong, ids, elms.n);
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

  slong *vrt = tcalloc(slong, elms.n * nv);
  for (uint i = 0; i < elms.n; i++) {
    for (int j = 0; j < nv; j++)
      vrt[i * nv + j] = pe[i].vtx[j];
  }

  if (verbose > 1) {
    for (uint i = 0; i < elms.n; i++) {
      printf("fid = %llu, ", ids[i]);
      for (int v = 0; v < nv; v++)
        printf("%lld, ", vrt[i * nv + v]);
      printf("\n");
      fflush(stdout);
    }
  }

  // Find and compress neighbors in order to form the Laplacian
  struct array nbrs, eij;
  find_nbrs(&nbrs, ids, vrt, elms.n, nv, cr, bfr);
  free(ids), free(vrt);
  compress_nbrs(&eij, &nbrs, bfr);
  array_free(&elms);
  array_free(&nbrs);

  // Setup the parallel CSR matrix
  par_csr_setup(&ilu->A, &eij, 0, bfr);
  array_free(&eij);

  if (verbose > 1) {
    for (int l = 0; l < ilu->nlvls; l++) {
      for (uint i = ilu->lvl_off[l]; i < ilu->lvl_off[l + 1]; i++) {
        printf("id = %d, lvl = %d e = %u\n", c->id, l + 1, ilu->A.rows[i]);
        fflush(stdout);
      }
    }
  }

  return 0;
}

struct ilu *ilu_setup(const uint n, const int nv, const long long *vtxll,
                      const int type, const double tol, const int iter,
                      const MPI_Comm comm, const int verbose) {
  struct comm c;
  comm_init(&c, comm);

  struct ilu *ilu = tcalloc(struct ilu, 1);
  crystal_init(&ilu->cr, &c);

  slong *vtx = tcalloc(slong, n * nv);
  for (uint i = 0; i < n * nv; i++)
    vtx[i] = vtxll[i];

  // Establish a numbering based on input
  slong out[2][1], buf[2][1], in = n;
  comm_scan(out, &c, gs_long, gs_add, &in, 1, buf);

  ulong *ids = tcalloc(ulong, n);
  for (uint i = 0; i < n; i++)
    ids[i] = out[0][0] + i + 1;

  buffer bfr;
  buffer_init(&bfr, 1024);

  uint *lvl_off = tcalloc(uint, 100 + n), *lvl_owner = lvl_off + 100;
  ulong *lvl_ids = tcalloc(ulong, n);
  int nlvls = find_lvls(lvl_off, lvl_owner, lvl_ids, n, nv, ids, vtx, 0,
                        &ilu->cr, verbose, &bfr);
  ilu_setup_aux(ilu, nlvls, lvl_off, lvl_owner, lvl_ids, n, nv, vtx, verbose,
                &bfr);

  char *val = getenv("PARRSB_DUMP_ILU_PRE");
  if (val != NULL && atoi(val) != 0)
    par_mat_dump("pre.txt", &ilu->A, &ilu->cr, &bfr);

  // Setup the ILU factors
  ilu->type = type, ilu->tol = tol, ilu->iter = iter;
  if (c.id == 0 && verbose > 0) {
    printf("iter = %d tol = %lf\n", ilu->iter, ilu->tol);
    fflush(stdout);
  }
  switch (type) {
  case 0:
    ilu0(ilu, &c, &bfr);
    break;
  case 1:
    iluti(ilu, &c, &bfr);
    break;
  default:
    break;
  }

  val = getenv("PARRSB_DUMP_ILU_POST");
  if (val != NULL && atoi(val) != 0)
    par_mat_dump("post.txt", &ilu->A, &ilu->cr, &bfr);

  free(ids), free(vtx), free(lvl_off), free(lvl_ids);
  buffer_free(&bfr), comm_free(&c);

  return ilu;
}

void ilu_free(struct ilu *ilu) {
  if (ilu) {
    crystal_free(&ilu->cr);
    if (ilu->nlvls > 0)
      par_mat_free(&ilu->A);
    if (ilu->lvl_off)
      free(ilu->lvl_off);
    free(ilu);
  }
}
