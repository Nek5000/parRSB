#include <stdlib.h>
#include <math.h>

#include <gencon-impl.h>
#include <sort.h>

static void initSegments(Mesh mesh, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  // Initialize globalId and ifSegment
  slong out[2][1], buf[2][1], in[1];
  in[0] = nPoints;
  comm_scan(out, c, gs_long, gs_add, in, 1, buf);
  slong start = out[0][0];
  uint i;
  for (i = 0; i < nPoints; i++) {
    points[i].ifSegment = 0;
    points[i].globalId = start + i;
  }

  // First rank with nPoints > 0 set ifSegment = 1
  sint rank = c->id;
  if (nPoints == 0)
    rank = c->np;

  comm_allreduce(c, gs_int, gs_min, &rank, 1, buf);

  if (c->id == rank)
    points[0].ifSegment = 1;
}

static int sortSegments(Mesh mesh, struct comm *c, int dim, buffer *bfr) {
  sint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  sint s = 0, e;
  while (s < nPoints) {
    // find the length of the segment
    for (e = s + 1; e < nPoints && points[e].ifSegment == 0; e++);

    // sort start to end based on dim
    switch (dim) {
      case 0:
        sarray_sort(struct Point_private, &points[s], e - s, x[0], 3, bfr);
        break;
      case 1:
        sarray_sort(struct Point_private, &points[s], e - s, x[1], 3, bfr);
        break;
      case 2:
        sarray_sort(struct Point_private, &points[s], e - s, x[2], 3, bfr);
        break;
      default:
        break;
    }

    sint i, sum = 0;
    for (i = s + 1; i < e; i++) {
        sum += points[i].ifSegment;
        points[i].ifSegment = 0;
    }

    if (sum > 0)
      points[s].ifSegment = 1;

    s = e;
  }

  return 0;
}

static int findLocalSegments(Mesh mesh, int i, GenmapScalar tolSquared) {
  Point pts = mesh->elements.ptr;
  sint npts = mesh->elements.n;
  int nDim = mesh->nDim;

  sint j;
  for (j = 1; j < npts; j++) {
    GenmapScalar d = sqrDiff(pts[j].x[i], pts[j - 1].x[i]);

    GenmapScalar dx = min(pts[j].dx, pts[j - 1].dx)*tolSquared;

    if (d > dx)
      pts[j].ifSegment = 1;
  }
}

static int mergeSegments(Mesh mesh, struct comm *c, int i, GenmapScalar tolSquared, int transfer, buffer *bfr) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  sint rank = c->id;
  sint size = c->np;

  struct Point_private lastp = points[nPoints-1];
  lastp.proc = (rank + 1)%size;

  struct array arr;
  array_init(struct Point_private, &arr, 1);
  array_cat(struct Point_private, &arr, &lastp, 1);

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct Point_private, &arr, proc, 1, &cr);

  uint n = arr.n;
  assert(n == 1);
  lastp = ((struct Point_private *) arr.ptr)[0];

  if (rank > 0) {
    GenmapScalar d = sqrDiff(lastp.x[i], points->x[i]);

    GenmapScalar dx = min(lastp.dx, points->dx)*tolSquared;

    if (d > dx)
      points->ifSegment = 1;
  }

  array_free(&arr);

  if (transfer > 0) {
    // If rank > 0, send i = 0,... n-1 where points[i].ifSegment == 0 to rank - 1
    n = 0;
    for (; n < nPoints && points[n].ifSegment == 0; n++)
      points[n].proc = rank - 1;
    for (; n < nPoints; n++)
      points[n].proc = rank;

    sarray_transfer(struct Point_private, &mesh->elements, proc, 0, &cr);
    crystal_free(&cr);

    sarray_sort(struct Point_private, mesh->elements.ptr, mesh->elements.n, globalId, 1, bfr);
  }

  return 0;
}

static slong countSegments(Mesh mesh, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  sint count = 0, i;
  for (i = 0; i < nPoints; i++)
    if (points[i].ifSegment > 0)
      count++;

  slong in, buf[2][1];
  in = count;
  comm_allreduce(c, gs_long, gs_add, &in, 1, buf);
  return in;
}

struct schedule {
  int dim;
  slong segments;
};

int findSegments(Mesh mesh, struct comm *c, GenmapScalar tol, int verbose, buffer *bfr) {
  // TODO: load balance
  parallel_sort(struct Point_private, &mesh->elements, x[0], genmap_gs_scalar, bin_sort, 0, c);

  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;
  int nDim = mesh->nDim;

  if (nDim == 3)
    sarray_sort_3(struct Point_private, points, nPoints, x[0], 3, x[1], 3, x[2], 3, bfr);
  else
    sarray_sort_2(struct Point_private, points, nPoints, x[0], 3, x[1], 3, bfr);

  comm_ext orig;
#ifdef MPI
  MPI_Comm_dup(c->c, &orig);
#endif

  struct comm nonZeroRanks;

  int bin = (nPoints > 0);
#ifdef MPI
  MPI_Comm new;
  MPI_Comm_split(orig, bin, c->id, &new);
  comm_init(&nonZeroRanks, new);
  MPI_Comm_free(&new);
#else
  comm_init(&nonZeroRanks, 1);
#endif

  struct schedule sched[3];
  GenmapScalar tolSquared = tol*tol;

  int t;
  for (t = 0; t < nDim; t++){
    sched[t].dim = t;
    initSegments(mesh, &nonZeroRanks);
    findLocalSegments(mesh, t, tolSquared);
    mergeSegments(mesh, &nonZeroRanks, t, tolSquared, 0, bfr);
    sched[t].segments = countSegments(mesh, &nonZeroRanks);
  }

  sarray_sort(struct schedule, sched, 3, segments, 1, bfr);

  initSegments(mesh, c);

  int d, merge = 1, err = 0;

  for (t = 0; t < nDim && err == 0; t++) {
    for (d = nDim - 1; d >= 0; d--) {
      int dim = sched[d].dim;
      if (bin > 0) {
        sortSegments(mesh, &nonZeroRanks, dim, bfr);
        findLocalSegments(mesh, dim, tolSquared);

        slong segments = countSegments(mesh, &nonZeroRanks);

        int rank = nonZeroRanks.id;
        if (rank == 0 && verbose > 0)
          printf("\tlocglob: %d %d %lld\n", t + 1, dim + 1, segments);

        mergeSegments(mesh, &nonZeroRanks, dim, tolSquared, merge, bfr);
        merge = 0;
      }

      comm_free(&nonZeroRanks);

      nPoints = mesh->elements.n;
      bin = (nPoints > 0);
#ifdef MPI
      MPI_Comm new;
      MPI_Comm_split(orig, bin, nonZeroRanks.id, &new);
      comm_init(&nonZeroRanks, new);
      MPI_Comm_free(&new);
#else
      comm_init(&nonZeroRanks, 1);
#endif
    }
  }

  comm_free(&nonZeroRanks);
#ifdef MPI
  MPI_Comm_free(&orig);
#endif

  return err;
}
