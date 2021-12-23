#include "genmap-impl.h"
#include "multigrid.h"

#define MIN(a, b) ((b) < (a) ? (b) : (a))

struct nbr_entry {
  ulong r, c;
  uint proc;
};

struct csr_entry {
  ulong r, c;
  GenmapScalar v;
};

struct vertex0 {
  GenmapULong elementId;
  GenmapULong vertexId;
  uint workProc;
};

static void find_neighbors(struct array *arr, struct rsb_element *elems,
                           sint nelt, int nv, struct comm *cc, buffer *buf) {
  slong out[2][1], bfr[2][1];
  slong lelt = nelt;
  comm_scan(out, cc, gs_long, gs_add, &lelt, 1, bfr);
  ulong elem_id = out[0][0] + 1;
  ulong sequence_id = elem_id * nv;

  struct array vertices;
  size_t size = nelt * nv;
  array_init(struct vertex0, &vertices, size);

  sint i, j;
  for (i = 0; i < lelt; i++) {
    for (j = 0; j < nv; j++) {
      struct vertex0 vrt = {.elementId = elem_id,
                            .vertexId = elems[i].vertices[j],
                            .workProc = elems[i].vertices[j] % cc->np};
      array_cat(struct vertex0, &vertices, &vrt, 1);
    }
    elem_id++;
  }

  struct crystal cr;
  crystal_init(&cr, cc);

  sarray_transfer(struct vertex0, &vertices, workProc, 1, &cr);
  size = vertices.n;
  struct vertex0 *vPtr = vertices.ptr;

  sarray_sort(struct vertex0, vPtr, size, vertexId, 1, buf);

  /* FIXME: Assumes quads or hexes */
  array_init(struct nbr_entry, arr, 10);
  sint s = 0, e;
  struct nbr_entry t;
  while (s < size) {
    e = s + 1;
    while (e < size && vPtr[s].vertexId == vPtr[e].vertexId)
      e++;
    int nnbrs = MIN(e, size) - s;

    for (i = s; i < MIN(e, size); i++) {
      t.r = vPtr[i].elementId;
      t.proc = vPtr[i].workProc;
      for (j = 0; j < nnbrs; j++) {
        t.c = vPtr[s + j].elementId;
        array_cat(struct nbr_entry, arr, &t, 1);
      }
    }
    s = e;
  }

  sarray_transfer(struct nbr_entry, arr, proc, 1, &cr);
  crystal_free(&cr);
  array_free(&vertices);
}

//------------------------------------------------------------------------------
// Laplacian - CSR
//
static struct gs_data *csr_gs_top(struct csr_laplacian *M, struct comm *c) {
  const uint rn = M->rn;
  const uint n = M->roff[rn];

  slong *ids;
  if (n > 0)
    GenmapMalloc(n, &ids);

  uint i, j;
  for (i = 0; i < rn; i++)
    for (j = M->roff[i]; j < M->roff[i + 1]; j++)
      if (M->rstart + i == M->col[j])
        ids[j] = M->col[j];
      else
        ids[j] = -M->col[j];

  struct gs_data *gsh = gs_setup(ids, n, c, 0, gs_pairwise, 0);

  if (n > 0)
    GenmapFree(ids);

  return gsh;
}

static void csr_init_aux(struct csr_laplacian *M, struct array *entries,
                         struct comm *c, buffer *buf) {
  uint i = 0;
  uint rn = 0;
  uint j;
  struct csr_entry *ptr = entries->ptr;
  while (i < entries->n) {
    j = i + 1;
    while (j < entries->n && ptr[i].r == ptr[j].r)
      j++;
    i = j;
    rn++;
  }

  M->rn = rn;
  if (M->rn == 0) {
    M->col = NULL;
    M->v = NULL;
    M->buf = NULL;
    M->diag = NULL;
    M->roff = NULL;
    return;
  } else {
    GenmapMalloc(entries->n, &M->col);
    GenmapMalloc(entries->n, &M->v);
    GenmapMalloc(entries->n, &M->buf);
    GenmapMalloc(M->rn, &M->diag);
    GenmapMalloc(M->rn + 1, &M->roff);
  }

  slong out[2][1], bf[2][1];
  slong in = M->rn;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bf);
  M->rstart = out[0][0] + 1;

  M->roff[0] = 0;
  rn = 0;
  i = 0;
  ptr = entries->ptr;
  while (i < entries->n) {
    M->col[i] = ptr[i].c;
    M->v[i] = ptr[i].v;

    double diag = 0.0;
    uint idx;
    if (ptr[i].r == ptr[i].c)
      idx = i;
    else
      diag = ptr[i].v;

    j = i + 1;
    while (j < entries->n && ptr[i].r == ptr[j].r) {
      M->col[j] = ptr[j].c;
      M->v[j] = ptr[j].v;

      if (ptr[j].r == ptr[j].c)
        idx = j;
      else
        diag += ptr[j].v;

      j++;
    }

    M->v[idx] = M->diag[rn++] = -diag;
    i = M->roff[rn] = j;
  }
  assert(rn == M->rn);

  M->gsh = csr_gs_top(M, c);
}

static int csr_init(struct laplacian *l, struct rsb_element *elems, uint lelt,
                    int nv, struct comm *c, buffer *buf) {
  struct array entries;
  find_neighbors(&entries, elems, lelt, nv, c, buf);

  struct nbr_entry *ptr = entries.ptr;
  sarray_sort_2(struct nbr_entry, ptr, entries.n, r, 1, c, 1, buf);

  struct array unique;
  array_init(struct csr_entry, &unique, entries.n);

  struct csr_entry t = {.r = ptr[0].r, .c = ptr[0].c, .v = -1.0};
  uint i;
  if (l->type & UNWEIGHTED) {
    for (i = 1; i < entries.n; i++)
      if (t.r != ptr[i].r || t.c != ptr[i].c) {
        array_cat(struct csr_entry, &unique, &t, 1);
        t.r = ptr[i].r;
        t.c = ptr[i].c;
      }
  } else if (l->type & WEIGHTED) {
    for (i = 1; i < entries.n; i++)
      if (t.r != ptr[i].r || t.c != ptr[i].c) {
        array_cat(struct csr_entry, &unique, &t, 1);
        t.r = ptr[i].r;
        t.c = ptr[i].c;
        t.v = -1.0;
      } else
        t.v -= 1.0;
  }
  array_cat(struct csr_entry, &unique, &t, 1);
  array_free(&entries);

  l->data = tcalloc(struct csr_laplacian, 1);
  csr_init_aux((struct csr_laplacian *)l->data, &unique, c, buf);

  array_free(&unique);

  return 0;
}

static int csr(GenmapScalar *v, struct laplacian *l, GenmapScalar *u,
               buffer *buf) {
  csr_mat_apply(v, (struct csr_laplacian *)l->data, u, buf);
  return 0;
}

static int csr_free(struct laplacian *l) {
  struct csr_laplacian *M = l->data;
  if (M != NULL)
    csr_mat_free(M);
  free(l->data);
  l->data = NULL;
}

//------------------------------------------------------------------------------
// Laplacian - GS
//

static int gs_weighted_init(struct laplacian *l, struct rsb_element *elems,
                            uint lelt, int nv, struct comm *c, buffer *buf) {

  uint npts = nv * lelt;
  slong *vertices = tcalloc(slong, npts);
  uint i, j;
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      vertices[i * nv + j] = elems[i].vertices[j];

  struct gs_laplacian *gl = l->data = tcalloc(struct gs_laplacian, 1);
  gl->u = tcalloc(GenmapScalar, npts);
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      gl->u[nv * i + j] = 1.0;

  gl->gsh = gs_setup(vertices, npts, c, 0, gs_crystal_router, 0);
  gs(gl->u, gs_double, gs_add, 0, gl->gsh, buf);

  gl->diag = tcalloc(GenmapScalar, lelt);
  for (i = 0; i < lelt; i++) {
    gl->diag[i] = 0.0;
    for (j = 0; j < nv; j++)
      gl->diag[i] += gl->u[nv * i + j];
  }

#if 0
  slong out[2][1], bf[2][1];
  slong in = lelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bf);
  slong row_start = out[0][0] + 1;

  int k;
  for (k = 0; k < c->np; k++) {
    genmap_barrier(c);
    if (c->id == k) {
      for (i = 0; i < lelt; i++)
        fprintf(stderr, "gs %lld: %.10lf\n", row_start + i, l->diag[i]);
    }
    fflush(stderr);
  }
#endif

  if (vertices != NULL)
    free(vertices);

  return 0;
}

static int gs_weighted(GenmapScalar *v, struct laplacian *l, GenmapScalar *u,
                       buffer *buf) {
  uint lelt = l->nel;
  int nv = l->nv;

  struct gs_laplacian *gl = l->data;

  GenmapInt i, j;
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      gl->u[nv * i + j] = u[i];

  gs(gl->u, gs_double, gs_add, 0, gl->gsh, buf);

  for (i = 0; i < lelt; i++) {
    v[i] = gl->diag[i] * u[i];
    for (j = 0; j < nv; j++)
      v[i] -= gl->u[nv * i + j];
  }

  return 0;
}

static int gs_weighted_free(struct laplacian *l) {
  struct gs_laplacian *gl = l->data;

  if (gl->u != NULL)
    free(gl->u);
  if (gl->diag != NULL)
    free(gl->diag);
  gs_free(gl->gsh);

  free(l->data);
  l->data = NULL;
}

//-----------------------------------------------------------------------------
// Laplacian - GPU
//
struct gpu_csr_entry {
  ulong r, c;
  GenmapScalar v;
  uint idx;
};

static int gpu_init_aux(struct gpu_laplacian *gl, struct array *unique,
                        uint lelt, int weighted) {
  gl->diag_ind = tcalloc(uint, lelt);
  gl->diag_val = tcalloc(GenmapScalar, lelt);
  gl->adj_off = tcalloc(uint, lelt);
  gl->adj_ind = tcalloc(uint, unique->n - lelt);
  gl->adj_val = NULL;
  if (weighted)
    gl->adj_val = tcalloc(GenmapScalar, unique->n - lelt);

  struct gpu_csr_entry *ptr = unique->ptr;
  uint rn = 0;
  GenmapScalar diag = 0.0;

  gl->adj_off[0] = 0;
  uint adjn = 0;
  if (ptr[0].r != ptr[0].c) {
    if (weighted)
      gl->adj_val[adjn] = ptr[0].v;
    gl->adj_ind[adjn++] = ptr[0].idx;
    diag += ptr[0].v;
  } else
    gl->diag_ind[rn] = ptr[0].idx;

  uint i;
  for (i = 1; i < unique->n; i++) {
    if (ptr[i - 1].r != ptr[i].r) {
      gl->diag_val[rn++] = diag, diag = 0.0;
      gl->adj_off[rn] = adjn;
    }
    if (ptr[i].r != ptr[i].c) {
      if (weighted)
        gl->adj_val[adjn] = ptr[i].v;
      gl->adj_ind[adjn++] = ptr[i].idx;
      diag += ptr[i].v;
    } else
      gl->diag_ind[rn] = ptr[i].idx;
  }
  gl->rn = rn;

  assert(rn == lelt); // Sanity check

  return 0;
}

static int gpu_init(struct laplacian *l, struct rsb_element *elems, uint lelt,
                    int nv, struct comm *c, buffer *buf) {
  struct array nbrs;
  find_neighbors(&nbrs, elems, lelt, nv, c, buf);
  sarray_sort(struct nbr_entry, nbrs.ptr, nbrs.n, c, 1, buf);
  assert(nbrs.n > 0); // Sanity check

  struct nbr_entry *ptr = nbrs.ptr;
  struct gpu_csr_entry t = {.r = ptr[0].r, .c = ptr[0].c, .v = -1.0, .idx = 0};

  struct array entries;
  array_init(struct gpu_csr_entry, &entries, 10);
  array_cat(struct gpu_csr_entry, &entries, &t, 1);

  struct gpu_laplacian *gl = l->data = tcalloc(struct gpu_laplacian, 1);

  uint max_col = lelt;
  gl->col_ids = tcalloc(ulong, max_col);
  gl->col_ids[0] = ptr[0].c;

  uint i, cid;
  for (i = 1, cid = 0; i < nbrs.n; i++) {
    t.c = ptr[i].c;
    t.r = ptr[i].r;
    if (ptr[i - 1].c != ptr[i].c) {
      t.idx = ++cid;
      if (cid == max_col) {
        max_col += max_col / 2 + 1;
        gl->col_ids = realloc(gl->col_ids, sizeof(ulong) * max_col);
      }
      gl->col_ids[cid] = ptr[i].c;
    }
    array_cat(struct gpu_csr_entry, &entries, &t, 1);
  }
  gl->cn = cid + 1;
  array_free(&nbrs);

  sarray_sort_2(struct gpu_csr_entry, entries.ptr, entries.n, r, 1, c, 1, buf);
  struct gpu_csr_entry *gptr = entries.ptr;
  t = gptr[0];

  // Compress the entries array
  struct array unique;
  array_init(struct gpu_csr_entry, &unique, 10);

  if (l->type & UNWEIGHTED) {
    for (i = 1; i < entries.n; i++)
      if (gptr[i].r != t.r || gptr[i].c != t.c) {
        array_cat(struct gpu_csr_entry, &unique, &t, 1);
        t = gptr[i];
      }
  } else if (l->type & WEIGHTED) {
    for (i = 1; i < entries.n; i++)
      if (gptr[i].r != t.r || gptr[i].c != t.c) {
        array_cat(struct gpu_csr_entry, &unique, &t, 1);
        t = gptr[i];
        t.v = -1.0;
      } else
        t.v -= 1.0;
  }
  array_cat(struct gpu_csr_entry, &unique, &t, 1);
  array_free(&entries);

  gpu_init_aux(gl, &unique, lelt, l->type & WEIGHTED);
  array_free(&unique);

  slong out[2][1], bf[2][1];
  slong in = lelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bf);
  slong row_start = out[0][0] + 1;

  slong *col_ids = tcalloc(slong, gl->cn);

  for (i = 0; gl->col_ids[i] < row_start; i++)
    col_ids[i] = -gl->col_ids[i];

  gl->ls = i;
  for (; i < gl->ls + lelt; i++)
    col_ids[i] = gl->col_ids[i];

  for (i = i + lelt; i < gl->cn; i++)
    col_ids[i] = -gl->col_ids[i];

  gl->gsh = gs_setup(col_ids, gl->cn, c, 0, gs_pairwise, 0);

  free(col_ids);

  return 0;
}

static int gpu_free(struct laplacian *l) {
  struct gpu_laplacian *gl = l->data;

  if (gl->adj_off != NULL)
    free(gl->adj_off);
  if (gl->adj_ind != NULL)
    free(gl->adj_ind);
  if (gl->adj_val != NULL)
    free(gl->adj_val);

  if (gl->diag_ind != NULL)
    free(gl->diag_ind);
  if (gl->diag_val != NULL)
    free(gl->diag_val);

  if (gl->col_ids != NULL)
    free(gl->col_ids);

  gs_free(gl->gsh);

  free(l->data);
  l->data = NULL;
}

//------------------------------------------------------------------------------
// Laplacian
//
int laplacian_init(struct laplacian *l, struct rsb_element *elems, uint nel,
                   int nv, int type, struct comm *c, buffer *buf) {
  l->type = type;
  l->nv = nv;
  l->nel = nel;

  if (type & CSR) {
    csr_init(l, elems, nel, nv, c, buf);
  } else if (type & GS) {
    assert(type & WEIGHTED);
    gs_weighted_init(l, elems, nel, nv, c, buf);
  } else if (type & GPU) {
    gpu_init(l, elems, nel, nv, c, buf);
  }

  return 0;
}

int laplacian(GenmapScalar *v, struct laplacian *l, GenmapScalar *u,
              buffer *buf) {
  if (l->type & CSR) {
    csr(v, l, u, buf);
  } else if (l->type & GS) {
    assert(l->type & WEIGHTED);
    gs_weighted(v, l, u, buf);
  } else if (l->type & GPU) {
    // gpu(v, l, u, buf);
  }

  return 0;
}

void laplacian_free(struct laplacian *l) {
  if (l->type & CSR)
    csr_free(l);
  else if (l->type & GS) {
    assert(l->type & WEIGHTED);
    gs_weighted_free(l);
  } else if (l->type & GPU)
    gpu_free(l);
}

#undef MIN
