#include <genmap-impl.h>
#include <genmap-multigrid.h>

#define MIN(a, b) ((b) < (a) ? (b) : (a))

struct nbr_entry {
  ulong r, c;
  uint proc;
};

struct csr_entry {
  ulong r, c;
  GenmapScalar v;
};

struct col {
  ulong c;
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
  array_init(vertex, &vertices, size);

  sint i, j;
  for (i = 0; i < lelt; i++) {
    for (j = 0; j < nv; j++) {
      vertex vrt = {.sequenceId = sequence_id,
                    .nNeighbors = 0,
                    .elementId = elem_id,
                    .vertexId = elems[i].vertices[j],
                    .workProc = elems[i].vertices[j] % cc->np};
      array_cat(vertex, &vertices, &vrt, 1);
      sequence_id++;
    }
    elem_id++;
  }

  struct crystal cr;
  crystal_init(&cr, cc);

  sarray_transfer(vertex, &vertices, workProc, 1, &cr);
  size = vertices.n;
  vertex *vPtr = vertices.ptr;

  sarray_sort(vertex, vPtr, size, vertexId, 1, buf);

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

#if 0
int laplacian_init(struct laplacian *l, struct rsb_element *elems, uint lelt,
                   int nv, struct comm *c, buffer *buf) {
  assert(lelt > 0);

  l->type = 0 | CSR | UNWEIGHTED;
  l->nel = lelt;
  l->nv = nv;

  struct array entries;
  find_neighbors(&entries, elems, lelt, nv, c, buf);

  slong out[2][1], bf[2][1];
  slong in = lelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bf);
  ulong start = out[0][0] + 1;

  struct array ucols;
  array_init(struct col, &ucols, 10);

  sarray_sort(struct csr_entry, entries.ptr, entries.n, c, 1, buf);
  struct csr_entry *ptr = entries.ptr;
  if (entries.n > 0) {
    array_cat(struct col, &ucols, &ptr[0], 1);
    ptr[0].idx = 0;
  }

  uint i;
  uint ncol = 0;
  for (i = 1; i < entries.n; i++) {
    if (ptr[i - 1].c != ptr[i].c) {
      array_cat(struct col, &ucols, &ptr[i].c, 1);
      ptr[i].idx = ++ncol;
    } else
      ptr[i].idx = ptr[i - 1].idx;
  }

  l->col_ids = tcalloc(ulong, ucols.n);
  memcpy(l->col_ids, ucols.ptr, sizeof(struct col) * ucols.n);

  sarray_sort_2(struct csr_entry, entries.ptr, entries.n, r, 1, c, 1, buf);
  ptr = entries.ptr;

  uint adjn = 0;
  for (i = 0; i < entries.n; i++)
    if (ptr[i].r != ptr[i].c)
      adjn++;

  l->adj_off = tcalloc(uint, lelt + 1);
  l->adj_ind = tcalloc(uint, adjn);
  adjn = 0;

  l->diag_ind = tcalloc(uint, lelt);
  l->diag_val = tcalloc(uint, lelt);
  uint diag_n = 0;

  uint rn = 0;
  l->adj_off[0] = 0;
  uint diag = 0;
  for (i = 0; i < entries.n; i++) {
    if (ptr[i].r != ptr[i].c)
      l->adj_ind[adjn++] = ptr[i].idx, diag++;
    else
      l->diag_ind[diag_n] = ptr[i].idx;
    if (ptr[i - 1].r != ptr[i].r) {
      l->adj_off[++rn] = adjn;
      l->diag[diag_n++] = diag;
    }
  }
  assert(rn == lelt);
  assert(diag_n == lelt);

  array_free(&ucols);
  array_free(&entries);

  return 0;
}
#endif

static struct gs_data *get_csr_top(struct csr_mat *M, struct comm *c) {
  const uint rn = M->rn;
  const uint n = M->roff[rn];

  slong *ids;
  if (n > 0)
    GenmapMalloc(n, &ids);

  uint i, j;
  for (i = 0; i < rn; i++)
    for (j = M->roff[i]; j < M->roff[i + 1]; j++)
      if (M->row_start + i == M->col[j])
        ids[j] = M->col[j];
      else
        ids[j] = -M->col[j];

  struct gs_data *gsh = gs_setup(ids, n, c, 0, gs_pairwise, 0);

  if (n > 0)
    GenmapFree(ids);

  return gsh;
}

static void csr_laplacian_setup(struct csr_mat *M, struct array *entries,
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
  M->row_start = out[0][0] + 1;

  M->roff[0] = 0;
  rn = 0;
  i = 0;
  ptr = entries->ptr;
  while (i < entries->n) {
    M->col[i] = ptr[i].c;
    M->v[i] = ptr[i].v;

    double diag = ptr[i].v;
    uint idx;
    if (ptr[i].r == ptr[i].c)
      idx = i;

    j = i + 1;
    while (j < entries->n && ptr[i].r == ptr[j].r) {
      M->col[j] = ptr[j].c;
      M->v[j] = ptr[j].v;

      diag += ptr[j].v;
      if (ptr[j].r == ptr[j].c)
        idx = j;

      j++;
    }

    M->v[idx] = M->diag[rn++] = -diag;
    i = M->roff[rn] = j;
  }
  assert(rn == M->rn);

  M->gsh = get_csr_top(M, c);
}

static int csr_unweighted_init(struct laplacian *l, struct rsb_element *elems,
                               uint lelt, int nv, struct comm *c, buffer *buf) {
  l->type = CSR | UNWEIGHTED;

  struct array entries;
  find_neighbors(&entries, elems, lelt, nv, c, buf);

  struct nbr_entry *ptr = entries.ptr;
  sarray_sort_2(struct nbr_entry, ptr, entries.n, r, 1, c, 1, buf);

  struct array unique;
  array_init(struct csr_entry, &unique, entries.n);

  struct csr_entry t = {.r = ptr[0].r, .c = ptr[0].c, .v = -1.0};
  array_cat(struct csr_entry, &unique, &t, 1);

  uint i;
  for (i = 1; i < entries.n; i++)
    if (t.r != ptr[i].r && t.c != ptr[i].c) {
      t.r = ptr[i].r;
      t.c = ptr[i].c;
      array_cat(struct csr_entry, &unique, &t, 1);
    }
  array_free(&entries);

  csr_laplacian_setup(l->M, &unique, c, buf);
  array_free(&unique);

  return 0;
}

static int csr_unweighted(GenmapScalar *v, struct laplacian *l, GenmapScalar *u,
                          buffer *buf) {
  csr_mat_apply(v, l->M, u, buf);
  return 0;
}

static int gs_weighted_init(struct laplacian *l, struct rsb_element *elems,
                            uint lelt, int nv, struct comm *c, buffer *buf) {
  l->type = GS | WEIGHTED;

  uint npts = nv * lelt;
  slong *vertices = tcalloc(slong, npts);
  uint i, j;
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      vertices[i * nv + j] = elems[i].vertices[j];

  l->u = tcalloc(GenmapScalar, npts);
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      l->u[nv * i + j] = 1.0;

  l->gsh = gs_setup(vertices, npts, c, 0, gs_crystal_router, 0);
  gs(l->u, gs_double, gs_add, 0, l->gsh, buf);

  l->diag = tcalloc(GenmapScalar, lelt);
  for (i = 0; i < lelt; i++) {
    l->diag[i] = 0.0;
    for (j = 0; j < nv; j++)
      l->diag[i] += l->u[nv * i + j];
  }

  if (vertices != NULL)
    free(vertices);

  return 0;
}

static int gs_weighted(GenmapScalar *v, struct laplacian *l, GenmapScalar *u,
                       buffer *buf) {
  uint lelt = l->nel;
  int nv = l->nv;

  GenmapInt i, j;
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      l->u[nv * i + j] = u[i];

  gs(l->u, gs_double, gs_add, 0, l->gsh, buf);

  for (i = 0; i < lelt; i++) {
    v[i] = l->diag[i] * u[i];
    for (j = 0; j < nv; j++)
      v[i] -= l->u[nv * i + j];
  }

  return 0;
}

int laplacian_init(struct laplacian *l, struct rsb_element *elems, uint nel,
                   int nv, int type, struct comm *c, buffer *buf) {
  l->type = type;

  l->diag = NULL;
  l->gsh = NULL;

  l->M = NULL;

  l->col_ids = NULL;
  l->adj_off = NULL;
  l->adj_ind = NULL;
  l->diag_ind = NULL;
  l->diag_val = NULL;

  l->u = NULL;

  if ((type & CSR) == CSR) {
    assert((type & UNWEIGHTED) == UNWEIGHTED);
    l->M = tcalloc(struct csr_mat, 1);
    csr_unweighted_init(l, elems, nel, nv, c, buf);
  } else if ((type & GS) == GS) {
    assert((type & WEIGHTED) == WEIGHTED);
    gs_weighted_init(l, elems, nel, nv, c, buf);
  }

  l->nv = nv;
  l->nel = nel;

  return 0;
}

int laplacian(GenmapScalar *v, struct laplacian *l, GenmapScalar *u,
              buffer *buf) {
  if ((l->type & CSR) == CSR) {
    assert((l->type & UNWEIGHTED) == UNWEIGHTED);
    csr_unweighted(v, l, u, buf);
  } else if ((l->type & GS) == GS) {
    assert((l->type & WEIGHTED) == WEIGHTED);
    gs_weighted(v, l, u, buf);
  }

  return 0;
}

void laplacian_free(struct laplacian *l) {
  if (l->diag != NULL)
    free(l->diag);
  if (l->gsh != NULL)
    gs_free(l->gsh);

  if (l->M != NULL)
    csr_mat_free(l->M);

  if (l->col_ids != NULL)
    free(l->col_ids);
  if (l->adj_off != NULL)
    free(l->adj_off);
  if (l->adj_ind != NULL)
    free(l->adj_ind);
  if (l->diag_ind != NULL)
    free(l->diag_ind);
  if (l->diag_val != NULL)
    free(l->diag_val);

  if (l->u != NULL)
    free(l->u);
}

#undef MIN
