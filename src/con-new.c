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

// Upper bounds for elements and face quantities
#define GC_MAX_FACES 6
#define GC_MAX_VERTICES 8
#define GC_MAX_NEIGHBORS 3
#define GC_MAX_FACE_VERTICES 4

// Varios mappings
static int PRE_TO_SYM_VERTEX[GC_MAX_VERTICES] = {0, 1, 3, 2, 4, 5, 7, 6};
static int PRE_TO_SYM_FACE[GC_MAX_FACES] = {2, 1, 3, 0, 4, 5};
static int NEIGHBOR_MAP[GC_MAX_VERTICES][GC_MAX_NEIGHBORS] = {
    {1, 2, 4}, {0, 3, 5}, {0, 3, 6}, {1, 2, 7},
    {0, 5, 6}, {1, 4, 7}, {2, 4, 7}, {3, 5, 6}};
static int SYM_FACES[GC_MAX_FACES][GC_MAX_FACE_VERTICES] = {
    {1, 3, 7, 5}, {2, 4, 8, 6}, // 1, 2
    {1, 2, 6, 5}, {3, 4, 8, 7}, // 3, 4
    {1, 2, 4, 3}, {5, 6, 8, 7}  // 5, 6
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

static void calc_centroid(struct face_t *face, unsigned nvf) {
  face->y[0] = face->y[1] = face->y[2] = 0;
  for (unsigned i = 0; i < nvf; i++) {
    face->y[0] += face->x[i][0];
    face->y[1] += face->x[i][1];
    face->y[2] += face->x[i][2];
  }
  face->y[0] /= nvf, face->y[1] /= nvf, face->y[2] /= nvf;
}

static void match_face_pair(struct face_t *f, struct face_t *g, unsigned nf,
                            unsigned nvf, unsigned nd) {
  scalar d2 = 0;
  for (unsigned i = 0; i < nvf; i++)
    d2 += D3d(f->x[i], g->x[i]);
  if (d2 < 1e-3 * f->dx && d2 < 1e-3 * g->dx) {
    f->gid[0] = g->gid[0] = MIN(f->gid[0], g->gid[0]);
    f->gid[1] = g->gid[1] = MIN(f->gid[1], g->gid[1]);
    f->gid[2] = g->gid[2] = MIN(f->gid[2], g->gid[2]);
    f->gid[3] = g->gid[3] = MIN(f->gid[3], g->gid[3]);
  }
}

static int match_faces(struct array *elems, ulong s, unsigned nv, unsigned nf,
                       unsigned nvf, unsigned ndim, struct crystal *cr,
                       buffer *bfr) {
  if (nv != 8 || ndim != 3 || nf != 6 || nvf != 4)
    return 1;

  struct comm *c = &cr->comm;

  struct array faces;
  array_init(struct face_t, &faces, nf * (elems->n / nv));

#define SET_FACE(pf, n, pe, e, f, v)                                           \
  {                                                                            \
    pf[n].dx = MIN(pf[n].dx, pe[e * nv + SYM_FACES[f][v] - 1].dx);             \
    pf[n].gid[v] = pe[e * nv + SYM_FACES[f][v] - 1].sid;                       \
    memcpy(pf[n].x[v], pe[e * nv + SYM_FACES[f][v] - 1].x,                     \
           sizeof(scalar) * ndim);                                             \
  }

  struct point_t *pe = (struct point_t *)elems->ptr;
  struct face_t *pf = (struct face_t *)faces.ptr;
  uint n = 0;
  for (uint e = 0; e < elems->n / nv; e++) {
    for (unsigned f = 0; f < nf; f++) {
      pf[n].dx = SCALAR_MAX;
      pf[n].p = c->id;
      pf[n].eid = pe[e * nv].eid;
      pf[n].sid = nf * pe[e * nv].eid + f;
      for (unsigned v = 0; v < nvf; v++)
        SET_FACE(pf, n, pe, e, f, v);
      calc_centroid(&pf[n], nvf);
      n++;
    }
  }
  assert(n == (elems->n / nv) * nf);
  faces.n = n;

#undef SET_FACE

  // Sort by centroid and match locally first
  sarray_sort_3(struct face_t, faces.ptr, faces.n, y[0], 3, y[1], 3, y[2], 3,
                bfr);

  if (faces.n > 0) {
    pf = (struct face_t *)faces.ptr;
    for (uint i = 1; i < faces.n; i++)
      match_face_pair(&pf[i - 1], &pf[i], nf, nvf, ndim);
  }

  sarray_sort_2(struct face_t, faces.ptr, faces.n, eid, 1, sid, 1, bfr);

  if (faces.n > 0) {
    pf = (struct face_t *)faces.ptr;
    for (uint n = 0; n < faces.n; n++) {
      ulong e = pf[n].eid - s;
      unsigned f = pf[n].sid % nf;
      for (unsigned v = 0; v < nvf; v++)
        pe[e * nv + SYM_FACES[f][v] - 1].sid =
            MIN(pe[e * nv + SYM_FACES[f][v] - 1].sid, pf[n].gid[v]);
    }
  }

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

  uint nelt = nelgt / c->np;
  uint nrem = nelgt - nelt * c->np;
  slong N = (c->np - nrem) * nelt;

  struct bc_t *pb = bcs->ptr;
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

  struct point_t *pe = (struct point_t *)elems.ptr;
  uint n = 0;
  for (uint i = 0; i < nelt; i++) {
    for (unsigned v = 0; v < nv; v++) {
      unsigned j = PRE_TO_SYM_VERTEX[v];
      for (unsigned k = 0; k < ndim; k++)
        pe[n].x[k] = coord[i * nv * ndim + j * ndim + k];
      pe[n].eid = s + i, pe[n].sid = nv * (s + i) + v, pe[n].p = c.id, n++;
    }
  }
  elems.n = n;

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
#undef GC_MAX_FACES
#undef GC_MAX_VERTICES
#undef GC_MAX_NEIGHBORS
#undef GC_MAX_FACE_VERTICES
