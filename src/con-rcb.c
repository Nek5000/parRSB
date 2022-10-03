#include "genmap-impl.h"
#include "sort.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Simple min/max macros
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// Upper bounds for elements and face quantities
#define GC_MAX_FACES 6
#define GC_MAX_VERTICES 8
#define GC_MAX_NEIGHBORS 3
#define GC_MAX_FACE_VERTICES 4

struct point_t {
  double dx, x[3];
  uint p, origin, ifseg;
  ulong seq_id, eid, gid;
};

struct periodic_t {
  ulong eid, fid;
  uint p;
  long long bc[2];
  struct point_t face[4];
};

//==============================================================================
// Transfer periodic faces in a load balance manner
//
static void transfer_bc_faces(struct array *pfaces, struct crystal *cr,
                              slong nelgt) {
  struct comm *c = &cr->comm;

  sint size = c->np;
  sint nelt = nelgt / size;
  sint nrem = nelgt - nelt * size;
  slong N = (size - nrem) * nelt;

  struct periodic_t *pf = (struct periodic_t *)pfaces->ptr;
  uint npf = pfaces->n;
  for (uint i = 0; i < npf; i++) {
    slong eid = pf[i].eid;
    if (eid < N)
      pf[i].p = eid / nelt;
    else
      pf[i].p = (eid - N) / (nelt + 1) + size - nrem;
  }

  sarray_transfer(struct periodic_t, pfaces, p, 1, cr);
}

//==============================================================================
// Distance functions
//
static inline double diff2(double x, double y) { return (x - y) * (x - y); }

static inline double distance2_2d(struct point_t *a, struct point_t *b) {
  return diff2(a->x[0], b->x[0]) + diff2(a->x[1], b->x[1]);
}

static inline double distance2_3d(struct point_t *a, struct point_t *b) {
  return distance2_2d(a, b) + diff2(a->x[2], b->x[2]);
}

//==============================================================================
// Find the minimum distance between a vertex and its neighbors
//
int NEIGHBOR_MAP[GC_MAX_VERTICES][GC_MAX_NEIGHBORS] = {
    {1, 2, 4}, {0, 3, 5}, {0, 3, 6}, {1, 2, 7},
    {0, 5, 6}, {1, 4, 7}, {2, 4, 7}, {3, 5, 6}};

static void find_dx(struct array *points, unsigned nd) {
  unsigned nv = (nd == 3 ? 8 : 4);
  unsigned nbrs = nd;

  struct point_t *p = (struct point_t *)points->ptr;
  if (nd == 3) {
    for (uint i = 0; i < points->n; i += nv) {
      for (unsigned j = 0; j < nv; j++) {
        p[i + j].dx = DBL_MAX;
        for (unsigned k = 0; k < nbrs; k++) {
          unsigned nbr = NEIGHBOR_MAP[j][k];
          double d = distance2_3d(&p[i + j], &p[i + nbr]);
          p[i + j].dx = MIN(p[i + j].dx, d);
        }
      }
    }
  } else if (nd == 2) {
    for (uint i = 0; i < points->n; i += nv) {
      for (unsigned j = 0; j < nv; j++) {
        p[i + j].dx = DBL_MAX;
        for (unsigned k = 0; k < nbrs; k++) {
          unsigned nbr = NEIGHBOR_MAP[j][k];
          double d = distance2_2d(&p[i + j], &p[i + nbr]);
          p[i + j].dx = MIN(p[i + j].dx, d);
        }
      }
    }
  }
}

//==============================================================================
// Perform RCB on vertex coordinates
//
static void bbox_local(double xmin[3], double xmax[3], struct array *points,
                       uint s, uint e, unsigned nd) {
  xmin[0] = xmin[1] = xmin[2] = DBL_MAX;
  xmax[0] = xmax[1] = xmax[2] = -DBL_MAX;

  struct point_t *p = (struct point_t *)points->ptr;
  if (nd == 3) {
    for (uint i = s; i < e; i++) {
      xmin[0] = (p[i].x[0] < xmin[0] ? p[i].x[0] : xmin[0]);
      xmin[1] = (p[i].x[1] < xmin[1] ? p[i].x[1] : xmin[1]);
      xmin[2] = (p[i].x[2] < xmin[2] ? p[i].x[2] : xmin[2]);
      xmax[0] = (p[i].x[0] > xmax[0] ? p[i].x[0] : xmax[0]);
      xmax[1] = (p[i].x[1] > xmax[1] ? p[i].x[1] : xmax[1]);
      xmax[2] = (p[i].x[2] > xmax[2] ? p[i].x[2] : xmax[2]);
    }
  } else if (nd == 2) {
    for (uint i = s; i < e; i++) {
      xmin[0] = (p[i].x[0] < xmin[0] ? p[i].x[0] : xmin[0]);
      xmin[1] = (p[i].x[1] < xmin[1] ? p[i].x[1] : xmin[1]);
      xmax[0] = (p[i].x[0] > xmax[0] ? p[i].x[0] : xmax[0]);
      xmax[1] = (p[i].x[1] > xmax[1] ? p[i].x[1] : xmax[1]);
    }
  }
}

static unsigned get_axis(const double xmin[3], const double xmax[3],
                         unsigned nd) {
  double wrk[3] = {xmax[0] - xmin[0], xmax[1] - xmin[1]};
  unsigned axis = (wrk[1] > wrk[0] ? 1 : 0);
  if (nd == 3) {
    wrk[2] = xmax[2] - xmin[2];
    axis = (wrk[2] > wrk[axis] ? 2 : axis);
  }

  return axis;
}

static void rcb_local(struct array *points, uint s, uint e, unsigned nd,
                      buffer *bfr) {
  if (s < e && e <= points->n) {
    double xmin[3], xmax[3];
    bbox_local(xmin, xmax, points, s, e, nd);
    unsigned axis = get_axis(xmin, xmax, nd);

    struct point_t *ps = (struct point_t *)points->ptr + s;
    uint size = e - s;
    switch (axis) {
    case 0:
      sarray_sort(struct point_t, ps, size, x[0], 3, bfr);
      break;
    case 1:
      sarray_sort(struct point_t, ps, size, x[1], 3, bfr);
      break;
    case 2:
      sarray_sort(struct point_t, ps, size, x[2], 3, bfr);
      break;
    default:
      break;
    }

    // TODO: May be replace this recursion with serial implementation
    rcb_local(points, s, (s + e) / 2, nd, bfr);
    rcb_local(points, (s + e) / 2, e, nd, bfr);
  }
}

static void rcb(struct array *points, unsigned nd, struct crystal *cr,
                buffer *bfr) {
  struct comm *ci = &cr->comm;
  struct comm c, t;

  comm_dup(&c, ci);
  while (c.np > 1) {
    // Find the axis for rcb
    double xmin[3], xmax[3], wrk[3];
    bbox_local(xmin, xmax, points, 0, points->n, nd);
    comm_allreduce(&c, gs_double, gs_min, xmin, 3, wrk);
    comm_allreduce(&c, gs_double, gs_max, xmax, 3, wrk);
    unsigned axis = get_axis(xmin, xmax, nd);

    // RCB on the identified axis
    switch (axis) {
    case 0:
      parallel_sort(struct point_t, points, x[0], gs_double, 0, 1, &c, bfr);
      break;
    case 1:
      parallel_sort(struct point_t, points, x[1], gs_double, 0, 1, &c, bfr);
      break;
    case 2:
      parallel_sort(struct point_t, points, x[2], gs_double, 0, 1, &c, bfr);
      break;
    }

    comm_split(&c, c.id < (c.np + 1) / 2, c.id, &t);
    comm_free(&c);
    comm_dup(&c, &t);
    comm_free(&t);
  }
  comm_free(&c);

  rcb_local(points, 0, points->n, nd, bfr);
}

//==============================================================================
// Resolve local connectivity
//
// Tuple sort
static void tuple_sort_(void *ra, uint n, uint usize, uint offset) {
  sint i, ir, j, l;
  void *rra = calloc(1, usize);
  assert(rra != NULL);

#define cpy(rra_, i_, ra_, l_)                                                 \
  {                                                                            \
    char *dst = (char *)rra_ + (i_ - 1) * usize;                               \
    char *src = (char *)ra_ + (l_ - 1) * usize;                                \
    memcpy(dst, src, usize);                                                   \
  }

#define get(ra_, l_) (*((double *)((char *)ra_ + (l_ - 1) * usize + offset)))

  if (n < 2)
    return;
  l = n / 2 + 1;
  ir = n;

  for (;;) {
    if (l > 1) {
      --l;
      assert(l >= 1 && l <= n && "l");
      cpy(rra, 1, ra, l);
    } else {
      cpy(rra, 1, ra, ir);
      cpy(ra, ir, ra, 1);
      if (--ir == 1) {
        cpy(ra, 1, rra, 1);
        break;
      }
    }
    i = l;
    j = l + l;
    while (j <= ir) {
      if (j < ir && get(ra, j) < get(ra, j + 1))
        j++;
      assert(j >= 1 && j <= n && "j2");
      assert(i >= 1 && i <= n && "i");
      if (get(rra, 1) < get(ra, j)) {
        cpy(ra, i, ra, j);
        i = j;
        j = 2 * j;
      } else
        break;
    }
    assert(i >= 1 && i <= n && "i2");
    cpy(ra, i, rra, 1);
  }

#undef cpy
#undef get

  if (rra != NULL)
    free(rra);
}

#define tuple_sort(T, arr, n, index)                                           \
  tuple_sort_((void *)arr, n, sizeof(T), offsetof(T, index))

static void local_connectivity(struct array *points, unsigned nd, double tol,
                               int verbose, buffer *bfr) {
  double tol2 = tol * tol;

  uint n = points->n;
  struct point_t *pp = (struct point_t *)points->ptr;
  if (n > 0) {
    // init segments
    pp[0].gid = 0, pp[0].ifseg = 1;
    for (uint i = 1; i < n; i++)
      pp[i].gid = pp[i].ifseg = 0;

    for (unsigned t = 0; t < nd; t++) {
      for (unsigned d = 0; d < nd; d++) {
        // sort segments
        uint s = 0, e;
        while (s < n - 1) {
          for (e = s + 1; e < n && pp[e].ifseg == 0; e++)
            ;

          if (e - s > 1) {
            switch (d) {
            case 0:
              tuple_sort(struct point_t, &pp[s], e - s, x[0]);
              break;
            case 1:
              tuple_sort(struct point_t, &pp[s], e - s, x[1]);
              break;
            case 2:
              tuple_sort(struct point_t, &pp[s], e - s, x[2]);
              break;
            default:
              break;
            }
          }

          s = e;
        }

        // find segments
        for (uint j = 1; j < n; j++) {
          double d2 = diff2(pp[j].x[d], pp[j - 1].x[d]);
          double dx2 = MIN(pp[j].dx, pp[j - 1].dx) * tol2;
          if (d2 > dx2)
            pp[j].ifseg = 1;
        }
      }
    }
  }

  // Now establish a local numbering of the points
  uint cnt = 0;
  for (uint i = 0; i < n; i++) {
    if (pp[i].ifseg)
      cnt++;
    pp[i].gid = cnt;
  }
}

//==============================================================================
// Set initial global numbering
//
static void initial_global_ids(struct array *points, struct comm *c) {
  uint n = points->n;
  slong out[2][1], wrk[2][1], in = n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong s = out[0][0];

  if (n > 0) {
    struct point_t *p = (struct point_t *)points->ptr;
    for (uint i = 0; i < n; i++)
      p[i].gid += s;
  }
}

//==============================================================================
// Resolve global ids going up the RCB tree
//
static void resolve_interface_ids(struct array *points, struct comm *c) {
  unsigned pcoord, nlvls;

  struct comm t;
  unsigned part = pcoord;
  for (unsigned i = 0; i < nlvls; i++) {
    // find the bounding box
    // find the interface nodes
    comm_dup(&t, c);
    // resolve interface nodes
    comm_free(&t);
    part /= 2;
  }
}

//==============================================================================
// C interface to find_conn
//
#define check_error(err, msg)                                                  \
  do {                                                                         \
    sint buf;                                                                  \
    comm_allreduce(&c, gs_int, gs_max, &err, 1, &buf);                         \
    if (err) {                                                                 \
      if (c.id == 0) {                                                         \
        printf("\n Error: %s\n", msg);                                         \
        fflush(stdout);                                                        \
      }                                                                        \
      buffer_free(&bfr);                                                       \
      comm_free(&c);                                                           \
      crystal_free(&cr);                                                       \
      return err;                                                              \
    }                                                                          \
  } while (0)

const int PRE_TO_SYM_VERTEX[GC_MAX_VERTICES] = {0, 1, 3, 2, 4, 5, 7, 6};
const int PRE_TO_SYM_FACE[GC_MAX_FACES] = {2, 1, 3, 0, 4, 5};

// Input:
//   nelt: Number of elements
//   ndim: Number of dimensions (2 or 3)
//   xyz [nelt, nv, ndim]: Coordinates of elements vertices in pre-processor
//     ordering, nv = 8 if ndim == 3 (Hex) or nv = 4 if ndim = 2 (Quad).
//   npf: Number of periodic faces
//   pf[npf, 4]: Periodic face information, 4 values: e1,f1,e2,f2
//   tol: Tolerance used to distinguish points
// Output:
//   vtx[nelt, nv]: Global numbering of vertices of elements
int parrsb_conn_rcb_mesh(long long *vtx, double *xyz, int nelt, int ndim,
                         long long *pf, int npf, double tol, MPI_Comm comm,
                         int verbose) {
  struct comm c;
  comm_init(&c, comm);

  struct crystal cr;
  crystal_init(&cr, &c);

  if (ndim != 2 && ndim != 3) {
    if (c.id == 0) {
      fprintf(stderr, "parCon: Dimesion of the mesh should be 2 or 3.\n");
      fflush(stderr);
    }
    return 1;
  }

  double duration[8] = {0};
  const char *name[8] = {"transfer_bc_faces     ", "find_dx               ",
                         "rcb                   ", "local_connectivity    ",
                         "initial_global_ids    ", "resolve_interface_ids ",
                         "faceCheck             ", "matchPeriodicFaces    "};

  if (c.id == 0 && verbose) {
    printf("Running parCon ... (tol = %g)\n", tol);
    fflush(stdout);
  }

  genmap_barrier(&c);
  double tcon = comm_time();

  slong out[2][1], wrk[2][1], in = nelt;
  comm_scan(out, &c, gs_long, gs_add, &in, 1, wrk);
  ulong start = out[0][0], nelgt = out[1][0];

  unsigned nv = (ndim == 3 ? 8 : 4);

  struct array points;
  array_init(struct point_t, &points, nelt * nv);

  struct point_t p = {.origin = c.id};
  for (uint i = 0; i < nelt; i++) {
    for (uint k = 0; k < nv; k++) {
      uint j = PRE_TO_SYM_VERTEX[k];
      for (uint l = 0; l < ndim; l++)
        p.x[l] = xyz[i * nv * ndim + j * ndim + l];
      p.eid = start + i;
      p.seq_id = nv * (start + i) + k;
      array_cat(struct point_t, &points, &p, 1);
    }
  }

  struct array pfaces;
  array_init(struct periodic_t, &pfaces, nelt * nv);

  struct periodic_t b;
  for (unsigned i = 0; i < npf; i++) {
    b.eid = pf[4 * i + 0] - 1;
    b.fid = PRE_TO_SYM_FACE[pf[4 * i + 1] - 1];
    b.bc[0] = pf[4 * i + 2] - 1;
    b.bc[1] = PRE_TO_SYM_FACE[pf[4 * i + 3] - 1];
    array_cat(struct periodic_t, &pfaces, &b, 1);
  }

  buffer bfr;
  buffer_init(&bfr, 1024);

  genmap_barrier(&c);
  double t = comm_time();
  transfer_bc_faces(&pfaces, &cr, nelgt);
  duration[0] = comm_time() - t;

  genmap_barrier(&c);
  t = comm_time();
  find_dx(&points, ndim);
  duration[1] = comm_time() - t;

  genmap_barrier(&c);
  t = comm_time();
  rcb(&points, ndim, &cr, &bfr);
  duration[2] = comm_time() - t;

  genmap_barrier(&c);
  t = comm_time();
  local_connectivity(&points, ndim, tol, verbose, &bfr);
  duration[3] = comm_time() - t;

  genmap_barrier(&c);
  t = comm_time();
  initial_global_ids(&points, &c);
  duration[4] = comm_time() - t;

  genmap_barrier(&c);
  t = comm_time();
  resolve_interface_ids(&points, &c);
  duration[5] = comm_time() - t;

#if 0
  genmap_barrier(&c);
  t = comm_time();
  err = elementCheck(mesh, &c, &bfr);
  check_error(err, "elementCheck");
  duration[4] = comm_time() - t;

  genmap_barrier(&c);
  t = comm_time();
  err = faceCheck(mesh, &c, &bfr);
  check_error(err, "faceCheck");
  duration[5] = comm_time() - t;

  genmap_barrier(&c);
  t = comm_time();
  err = matchPeriodicFaces(mesh, &c, &bfr);
  check_error(err, "matchPeriodicFaces");
  duration[6] = comm_time() - t;

  // Copy output
  Point ptr = mesh->points.ptr;
  for (i = 0; i < nelt; i++) {
    for (j = 0; j < nv; j++)
      vtx[i * nv + j] = ptr[i * nv + j].globalId + 1;
  }

  // Report timing info and finish
  genmap_barrier(&c);
  tcon = comm_time() - tcon;
  if (c.id == 0 && verbose > 0) {
    printf("parCon finished in %g s\n", tcon);
    fflush(stdout);
  }

  double gmin[8], gmax[8], buf[8];
  for (unsigned i = 0; i < 8; i++)
    gmax[i] = gmin[i] = duration[i];
  comm_allreduce(&c, gs_double, gs_min, gmin, 8, buf);
  comm_allreduce(&c, gs_double, gs_max, gmax, 8, buf);

  if (c.id == 0 && verbose > 0) {
    for (unsigned i = 0; i < 7; i++)
      printf("%s: %e %e (min max)\n", name[i], gmin[i], gmax[i]);
    fflush(stdout);
  }
#endif

  array_free(&points);
  array_free(&pfaces);
  buffer_free(&bfr);
  crystal_free(&cr);
  comm_free(&c);

  return 0;
}

//=============================================================================
// Fortran interface
//
void fparrsb_conn_rcb_mesh(long long *vtx, double *xyz, int *nelt, int *ndim,
                           long long *pf, int *npf, double *tol,
                           MPI_Fint *fcomm, int *verbose, int *err) {
  *err = 1;
  MPI_Comm c = MPI_Comm_f2c(*fcomm);
  *err =
      parrsb_conn_rcb_mesh(vtx, xyz, *nelt, *ndim, pf, *npf, *tol, c, *verbose);
}

#undef GC_MAX_FACES
#undef GC_MAX_VERTICES
#undef GC_MAX_NEIGHBORS
#undef GC_MAX_FACE_VERTICES
#undef MIN
#undef MAX
