#include "gslib.h"
#include <float.h>

#ifdef scalar
#undef scalar
#endif
#define scalar double

#ifdef SCALAR_MAX
#undef SCALAR_MAX
#endif
#define SCALAR_MAX DBL_MAX

#ifdef gs_scalar
#undef gs_scalar
#endif
#define gs_scalar gs_double

#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define CHK_ERR(err, msg, c)                                                   \
  {                                                                            \
    sint wrk;                                                                  \
    comm_allreduce((c), gs_int, gs_max, &err, 1, &wrk);                        \
    if (err != 0) {                                                            \
      if ((c)->id == 0)                                                        \
        printf("Error: %s:%d %s\n", __FILE__, __LINE__, msg);                  \
      return err;                                                              \
    }                                                                          \
  }

// Various mappings
static int PRE_TO_SYM_VERTEX[8] = {0, 1, 3, 2, 4, 5, 7, 6};

static int PRE_TO_SYM_FACE[6] = {2, 1, 3, 0, 4, 5};

static int NEIGHBOR_MAP[8][3] = {{1, 2, 4}, {0, 3, 5}, {0, 3, 6}, {1, 2, 7},
                                 {0, 5, 6}, {1, 4, 7}, {2, 4, 7}, {3, 5, 6}};

static int SYM_FACES[6][4] = {
    {0, 2, 4, 6}, {1, 3, 5, 7}, // 1, 2
    {0, 1, 4, 5}, {2, 3, 6, 7}, // 3, 4
    {0, 1, 2, 3}, {4, 5, 6, 7}  // 5, 6
};

struct point_t {
  scalar dx, x[3];
  uint p;
  ulong eid, sid;
};

static inline double D2(double x, double y) { return (x - y) * (x - y); }

static inline double D2d(double *x, double *y) {
  return D2(x[0], y[0]) + D2(x[1], y[1]);
}

static inline double D3d(double x[3], double y[3]) {
  return D2(x[0], y[0]) + D2(x[1], y[1]) + D2(x[2], y[2]);
}

struct bc_t {
  ulong bc[4];
  uint proc;
};

struct face_t {
  ulong eid, sid, gid[4];
  uint p;
  scalar x[4][3], y[3], dx;
};

// FIXME: Fix the calculation of the centroids: Use the mid point connecting
// the centroid of the triangles.
static void calc_centroid(struct face_t *face, unsigned nvf) {
  face->y[0] = face->y[1] = face->y[2] = 0;
  for (unsigned i = 0; i < nvf; i++) {
    face->y[0] += face->x[i][0];
    face->y[1] += face->x[i][1];
    face->y[2] += face->x[i][2];
  }
  face->y[0] /= nvf, face->y[1] /= nvf, face->y[2] /= nvf;
}

static void sorted_idx(unsigned idx[4], struct face_t *f, unsigned nvf,
                       unsigned nd, buffer *bfr) {
  struct idx_t {
    scalar x[3];
    uint idx;
  };

  struct array fpts;
  array_init(struct idx_t, &fpts, nvf);

  struct idx_t ti;
  for (unsigned i = 0; i < nvf; i++) {
    memcpy(ti.x, f->x[i], sizeof(scalar) * nd);
    ti.idx = i;
    array_cat(struct idx_t, &fpts, &ti, 1);
  }

  sarray_sort_3(struct idx_t, fpts.ptr, fpts.n, x[0], 3, x[1], 3, x[2], 3, bfr);
  if (fpts.n > 0) {
    struct idx_t *pf = (struct idx_t *)fpts.ptr;
    for (unsigned i = 0; i < fpts.n; i++)
      idx[i] = pf[i].idx;
  }

  array_free(&fpts);
}

static void match_face_pair(struct face_t *f, struct face_t *g, unsigned nvf,
                            unsigned nd, buffer *bfr) {
  // Find the sorted indices of the points of a face
  unsigned fidx[4] = {0}, gidx[4] = {0};
  sorted_idx(fidx, f, nvf, nd, bfr);
  sorted_idx(gidx, g, nvf, nd, bfr);

  // Find the pointwise distance now
  scalar d2 = 0;
  for (unsigned i = 0; i < nvf; i++)
    d2 += D3d(f->x[fidx[i]], g->x[gidx[i]]);

  if (d2 < 1e-3 * f->dx && d2 < 1e-3 * g->dx) {
    // It's a match
    for (unsigned i = 0; i < nvf; i++) {
      ulong min = MIN(f->gid[fidx[i]], g->gid[gidx[i]]);
      f->gid[fidx[i]] = g->gid[gidx[i]] = min;
    }
  }
}

static int match_faces(struct array *points, ulong s, unsigned nv, unsigned nf,
                       unsigned nvf, unsigned nd, struct crystal *cr,
                       buffer *bfr) {
  if (nv != 8 || nd != 3 || nf != 6 || nvf != 4)
    return 1;

  struct array faces;
  array_init(struct face_t, &faces, nf * (points->n / nv));

#define SET_FACE(fb, pfv)                                                      \
  {                                                                            \
    fb.dx = MIN(fb.dx, pfv->dx);                                               \
    fb.gid[v] = pfv->sid;                                                      \
    memcpy(fb.x[v], pfv->x, sizeof(scalar) * nd);                              \
  }

  struct comm *c = &cr->comm;
  struct point_t *pp = (struct point_t *)points->ptr;
  struct face_t fb = {.dx = SCALAR_MAX, .p = c->id};
  for (uint e = 0; e < points->n / nv; e++) {
    fb.eid = pp[e * nv].eid;
    for (unsigned f = 0; f < nf; f++) {
      fb.sid = nf * fb.eid + f;
      for (unsigned v = 0; v < nvf; v++) {
        struct point_t *pfv = &pp[e * nv + SYM_FACES[f][v]];
        SET_FACE(fb, pfv);
      }
      calc_centroid(&fb, nvf);
      array_cat(struct face_t, &faces, &fb, 1);
    }
  }

#undef SET_FACE

  ulong *initial = tcalloc(ulong, points->n);
  ulong *updated = tcalloc(ulong, points->n);
  for (uint i = 0; i < points->n; i++)
    initial[i] = updated[i] = pp[i].sid;

  int changed;
  do {
    // Sort by centroid and match locally first
    sarray_sort_3(struct face_t, faces.ptr, faces.n, y[0], 3, y[1], 3, y[2], 3,
                  bfr);

    if (faces.n > 0) {
      struct face_t *pf = (struct face_t *)faces.ptr;
      for (uint i = 1; i < faces.n; i++)
        match_face_pair(&pf[i - 1], &pf[i], nvf, nd, bfr);
    }

    sarray_sort_2(struct face_t, faces.ptr, faces.n, eid, 1, sid, 1, bfr);

    if (faces.n > 0) {
      struct face_t *pf = (struct face_t *)faces.ptr;
      for (uint n = 0; n < faces.n; n++) {
        ulong e = pf[n].eid - s;
        unsigned f = pf[n].sid % nf;
        for (unsigned v = 0; v < nvf; v++) {
          struct point_t *pfv = &pp[e * nv + SYM_FACES[f][v]];
          ulong min = MIN(pfv->sid, pf[n].gid[v]);
          pfv->sid = min;
        }
      }
    }

    changed = 0;
    for (uint i = 0; i < points->n; i++) {
      changed += (updated[i] != pp[i].sid);
      updated[i] = pp[i].sid;
    }
  } while (changed);

  if (points->n > 0) {
    for (uint i = 0; i < points->n; i++)
      printf("AFTER eid = %lld sid = %lld\n", pp[i].eid, pp[i].sid);
  }

  // If, updated == initialized, the connectivity never changed. This could be a
  // cornor vertex or two elements only connected by vertex or an edge without
  // any neighnors (or in other words, no other path to reach each other without
  // that edge or vertex).

  if (initial)
    free(initial);
  if (updated)
    free(updated);
  array_free(&faces);

  return 0;
}

static int find_min_nbr_distance(struct array *elems, unsigned ndim,
                                 unsigned nbrs, unsigned nv) {
  if (elems->n > 0) {
    struct point_t *pe = (struct point_t *)elems->ptr;
    if (ndim == 3) {
      for (uint i = 0; i < elems->n / nv; i++) {
        for (unsigned j = 0; j < nv; j++) {
          struct point_t *p = &pe[i * nv + j];
          p->dx = SCALAR_MAX;
          for (unsigned k = 0; k < nbrs; k++) {
            unsigned nbr = NEIGHBOR_MAP[j][k];
            p->dx = MIN(p->dx, D3d(p->x, pe[i * nv + nbr].x));
          }
        }
      }
    } else if (ndim == 2) {
      for (uint i = 0; i < elems->n / nv; i++) {
        for (unsigned j = 0; j < nv; j++) {
          struct point_t *p = &pe[i * nv + j];
          p->dx = SCALAR_MAX;
          for (unsigned k = 0; k < nbrs; k++) {
            unsigned nbr = NEIGHBOR_MAP[j][k];
            p->dx = MIN(p->dx, D2d(p->x, pe[i * nv + nbr].x));
          }
        }
      }
    }
  }

  return 0;
}

static int transfer_bc_faces(struct array *bcs, ulong nelgt,
                             struct crystal *cr) {
  struct comm *c = &cr->comm;
  uint nelt = nelgt / c->np, nrem = nelgt - nelt * c->np;

  slong N = (c->np - nrem) * nelt;
  struct bc_t *pb = (struct bc_t *)bcs->ptr;
  for (uint i = 0; i < bcs->n; i++) {
    slong eid = pb[i].bc[0];
    if (eid < N)
      pb[i].proc = eid / nelt;
    else
      pb[i].proc = (eid - N) / (nelt + 1) + c->np - nrem;
  }

  sarray_transfer(struct bc_t, bcs, proc, 1, cr);

  return 0;
}

// Input:
//   nelt: Number of elements, nv: Number of vertices in an element
//   coord [nelt, nv, ndim]: Coordinates of elements vertices in preprocessor
//     ordering, nv = 8 if ndim == 3 (Hex) or nv = 4 if ndim = 2 (Quad).
// Output:
//   vtx[nelt, nv]: Global numbering of vertices of elements
int parrsb_conn_new(long long *vtx, double *coord, int nelt, int ndim,
                    long long *pinfo, int npfaces, double tol, MPI_Comm comm,
                    int verbose) {
  struct comm c;
  comm_init(&c, comm);

  struct crystal cr;
  crystal_init(&cr, &c);

  if (c.id == 0 && verbose) {
    printf("Running parCon ... (tol=%g)\n", tol);
    fflush(stdout);
  }

  slong out[2][1], wrk[2][1], in = nelt;
  comm_scan(out, &c, gs_long, gs_add, &in, 1, wrk);
  ulong s = out[0][0], ng = out[1][0];

  unsigned nv = (ndim == 3 ? 8 : 4);
  struct array elems;
  array_init(struct point_t, &elems, nelt * nv);

  struct point_t pnt = {.p = c.id};
  for (uint i = 0; i < nelt; i++) {
    pnt.eid = s + i;
    for (unsigned v = 0; v < nv; v++) {
      unsigned j = PRE_TO_SYM_VERTEX[v];
      for (unsigned k = 0; k < ndim; k++)
        pnt.x[k] = coord[i * nv * ndim + j * ndim + k];
      pnt.sid = nv * pnt.eid + v;
      array_cat(struct point_t, &elems, &pnt, 1);
    }
  }

  struct array bcs;
  array_init(struct bc_t, &bcs, npfaces);
  struct bc_t *pb = (struct bc_t *)bcs.ptr;
  for (uint i = 0; i < npfaces; i++) {
    pb[i].bc[0] = pinfo[4 * i + 0] - 1;
    pb[i].bc[1] = PRE_TO_SYM_FACE[pinfo[4 * i + 1] - 1];
    pb[i].bc[2] = pinfo[4 * i + 2] - 1;
    pb[i].bc[3] = PRE_TO_SYM_FACE[pinfo[4 * i + 3] - 1];
  }
  bcs.n = npfaces;

  if (elems.n > 0) {
    struct point_t *pe = (struct point_t *)elems.ptr;
    for (uint i = 0; i < elems.n; i++)
      printf("BEFORE eid = %lld sid = %lld\n", pe[i].eid, pe[i].sid);
  }

  int err = transfer_bc_faces(&bcs, ng, &cr);
  CHK_ERR(err, "transfer_bc_faces", &c);

  err = find_min_nbr_distance(&elems, ndim, ndim, nv);
  CHK_ERR(err, "find_min_nbr_distance", &c);

  buffer bfr;
  buffer_init(&bfr, 1024);

  err = match_faces(&elems, s, nv, 2 * ndim, nv / 2, ndim, &cr, &bfr);
  CHK_ERR(err, "match_faces", &c);

  buffer_free(&bfr), array_free(&bcs), array_free(&elems);
  comm_free(&c), crystal_free(&cr);

  return 0;
}

#undef scalar
#undef SCALAR_MAX
#undef gs_scalar
#undef MIN
#undef CHK_ERR
