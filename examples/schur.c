//=============================================================================
// Test Schur complement solver
//
#include "coarse.h"
#include "parRSB.h"

#include <math.h>
#include <time.h>

static void schur_test(const unsigned int nelt, const int nv,
                       const long long *vl, MPI_Comm comm) {
  int id, np;
  MPI_Comm_rank(comm, &id);
  MPI_Comm_size(comm, &np);

  if (id == 0) {
    printf("MPI Ranks = %d", np);
    fflush(stdout);
  }

  MPI_Barrier(comm);
  double t = MPI_Wtime();
  struct coarse *crs = coarse_setup(nelt, nv, vl, 0, comm);
  t = MPI_Wtime() - t;

  if (id == 0) {
    printf("coarse_setup: %lf\n", t);
    fflush(stdout);
  }

  srand(time(NULL));

  double *b = (double *)tcalloc(double, nelt);
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
    b[i] /= sum;

  double *x = tcalloc(double, nelt);

  buffer bfr;
  buffer_init(&bfr, 1024);

  MPI_Barrier(comm);
  t = MPI_Wtime();
  coarse_solve(x, b, crs, &bfr);
  t = MPI_Wtime() - t;

  if (id == 0) {
    printf("coarse_solve: %lf\n", t);
    fflush(stdout);
  }

  free(b), free(x);
  buffer_free(&bfr);
  coarse_free(crs);
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

  schur_test(nelt, nv, vl, comm);

  // Free resources
  free(vl), free(coord);
  free(in);
  MPI_Finalize();

  return 0;
}
