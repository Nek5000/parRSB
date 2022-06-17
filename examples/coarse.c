//=============================================================================
// Test Schur complement solver
//
#include "coarse.h"
#include "parRSB.h"

#include <math.h>
#include <time.h>

static double check_err(double *b, double *x, uint nelt, int nv, slong *vtx,
                        MPI_Comm comm, buffer *bfr) {
  struct comm c;
  comm_init(&c, comm);

  struct crystal cr;
  crystal_init(&cr, &c);

  slong out[2][1], buf[2][1], in = nelt;
  comm_scan(out, &c, gs_long, gs_add, &in, 1, buf);
  ulong start = out[0][0] + 1;

  ulong *eid = tcalloc(ulong, nelt);
  for (uint i = 0; i < nelt; i++)
    eid[i] = start + i;

  struct array nbrs, eij;
  find_nbrs(&nbrs, eid, vtx, nelt, nv, &cr, bfr);
  compress_nbrs(&eij, &nbrs, bfr);

  struct par_mat M;
  par_csr_setup(&M, &eij, 1, bfr);
  assert(M.rn > 0);

  free(eid), array_free(&nbrs), array_free(&eij);

  struct gs_data *gsh = setup_Q(&M, &c, bfr);
  double *bl = tcalloc(double, nelt);
  double *wrk = tcalloc(double, M.rn + M.adj_off[M.rn]);
  mat_vec_csr(bl, x, &M, gsh, wrk, bfr);

  crystal_free(&cr), comm_free(&c);
  gs_free(gsh), par_mat_free(&M);

  double norm = 0.0;
  for (uint i = 0; i < nelt; i++)
    norm += (bl[i] - b[i]) * (bl[i] - b[i]);
  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, comm);

  free(wrk), free(bl);

  return sqrt(norm);
}

static void setup_rhs(double *b, const unsigned int nelt, MPI_Comm comm) {
  srand(time(NULL));
  double sum = 0;
  for (int i = 0; i < nelt; i++) {
    b[i] = (rand() % 50 + 1.0) / 10;
    sum += b[i];
  }
  MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);

  long long ng = nelt;
  MPI_Allreduce(MPI_IN_PLACE, &ng, 1, MPI_LONG_LONG, MPI_SUM, comm);
  sum /= ng;

  double norm = 0;
  for (int i = 0; i < nelt; i++) {
    b[i] -= sum;
    norm += b[i] * b[i];
  }

  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, comm);
  norm = sqrt(norm);

  for (int i = 0; i < nelt; i++)
    b[i] /= norm;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  struct parrsb_input *in = parrsb_parse_input(argc, argv, comm);
  int err = (in == NULL);
  parrsb_check_error(err, comm);

  // Read the geometry from the .re2 file, find connectiviy, partition and then
  // distribute the mesh.
  unsigned int nelt;
  int nv;
  long long *vl = NULL;
  double *coord = NULL;
  parrsb_setup_mesh(&nelt, &nv, &vl, &coord, in, comm);

  int ndim = (nv == 8) ? 3 : 2;
  double *centroids = tcalloc(double, nelt *ndim);
  for (uint i = 0; i < nelt; i++) {
    for (int j = 0; j < nv; j++) {
      for (int d = 0; d < ndim; d++)
        centroids[i * ndim + d] += coord[i * ndim * nv + j * ndim + d];
    }
    for (int d = 0; d < ndim; d++)
      centroids[i * ndim + d] /= nv;
  }

  // Setup the coarse solve with schur complement solver
  int id, np;
  MPI_Comm_rank(comm, &id);
  MPI_Comm_size(comm, &np);

  MPI_Barrier(comm);
  double t = MPI_Wtime();
  struct coarse *crs =
      coarse_setup(nelt, nv, vl, centroids, in->crs_type, comm);
  double tsetup = MPI_Wtime() - t;

  double *b = tcalloc(double, nelt);
  double *x = tcalloc(double, nelt);
  setup_rhs(b, nelt, comm);

  buffer bfr;
  buffer_init(&bfr, 1024);
  MPI_Barrier(comm);
  t = MPI_Wtime();
  coarse_solve(x, b, in->crs_tol, crs, &bfr);
  double tsolve = MPI_Wtime() - t;

  double enorm = check_err(b, x, nelt, nv, vl, comm, &bfr);
  if (id == 0) {
    printf("MPI Ranks = %d\ncoarse_setup: %lf\ncoarse_solve = %lf\nerr = %lf\n",
           np, tsetup, tsolve, enorm);
    fflush(stdout);
  }
  err = (enorm > 10 * in->crs_tol);
  parrsb_check_error(err, comm);

  // Free resources
  buffer_free(&bfr);
  coarse_free(crs);
  free(b), free(x), free(vl), free(coord), free(centroids);
  free(in);
  MPI_Finalize();

  return 0;
}
