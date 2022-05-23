//=============================================================================
// Test Schur complement solver
//
#include "coarse.h"
#include "parRSB.h"

#include <math.h>
#include <time.h>

static void schur_test(const unsigned int nelt, const int nv,
                       const long long *vl, MPI_Comm comm) {
  int id;
  MPI_Comm_rank(comm, &id);

  MPI_Barrier(comm);
  double t = MPI_Wtime();
  struct coarse *crs = coarse_setup(nelt, nv, vl, comm);
  t = MPI_Wtime() - t;

  if (id == 0) {
    printf("coarse_setup: %lf\n", t);
    fflush(stdout);
  }

  srand(time(NULL));

  double *b = (double *)tcalloc(double, nelt);
  double sum = 0;
  for (int i = 0; i < nelt; i++) {
    b[i] = (rand() % 50 + 1.0) / 2.5;
    sum += b[i];
  }
  MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);

  long ng = nelt;
  MPI_Allreduce(MPI_IN_PLACE, &ng, 1, MPI_LONG, MPI_SUM, comm);

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

  free(b), free(x), buffer_free(&bfr), coarse_free(crs);
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  struct parrsb_input *in = parrsb_parse_input(argc, argv, MPI_COMM_WORLD);
  int err = (in == NULL);
  parrsb_check_error(err, MPI_COMM_WORLD);

  int id;
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  // Read the geometry from the .re2 file
  unsigned int nelt, nbcs;
  double *coord = NULL;
  long long *bcs = NULL;
  int nv;
  err = parrsb_read_mesh(&nelt, &nv, NULL, &coord, &nbcs, &bcs, in->mesh,
                         MPI_COMM_WORLD, 1);
  parrsb_check_error(err, MPI_COMM_WORLD);

  // Find connectivity
  long long *vl = (long long *)calloc(nelt * nv, sizeof(long long));
  err = (vl == NULL);
  parrsb_check_error(err, MPI_COMM_WORLD);

  int ndim = (nv == 8 ? 3 : 2);
  err = parrsb_find_conn(vl, coord, nelt, ndim, bcs, nbcs, in->tol,
                         MPI_COMM_WORLD, 0);
  parrsb_check_error(err, MPI_COMM_WORLD);

  // Partition the mesh
  int *part = (int *)calloc(nelt, sizeof(int));
  err = (part == NULL);
  parrsb_check_error(err, MPI_COMM_WORLD);

  parrsb_options options = parrsb_default_options;
  err = parrsb_part_mesh(part, NULL, vl, coord, nelt, nv, options,
                         MPI_COMM_WORLD);
  parrsb_check_error(err, MPI_COMM_WORLD);

  // Redistribute data based on identified partitions
  err =
      parrsb_distribute_elements(&nelt, &vl, &coord, part, nv, MPI_COMM_WORLD);
  parrsb_check_error(err, MPI_COMM_WORLD);

  schur_test(nelt, nv, vl, MPI_COMM_WORLD);

  // Free resources
  free(part), free(vl), free(coord), free(bcs), free(in);
  MPI_Finalize();

  return 0;
}
