/*
 * Generate connectivity (.co2) from Nek5000 mesh file (.re2).
 */
#include <mpi.h>
#include <stdio.h>

#include <parRSB.h>

void check_error_(int err, char *file, int line, MPI_Comm comm) {
  int sum;
  MPI_Allreduce(&err, &sum, 1, MPI_INT, MPI_SUM, comm);

  if (sum != 0) {
    int id;
    MPI_Comm_rank(comm, &id);
    if (id == 0)
      printf("check_error failure in %s:%d\n", file, line);

    MPI_Finalize();
    exit(1);
  }
}

#define check_error(err) check_error_(err, __FILE__, __LINE__, MPI_COMM_WORLD);

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int id;
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  if (argc < 2 || argc > 3) {
    if (id == 0)
      printf("Usage: ./%s <mesh file> [tol]\n", argv[0]);
    MPI_Finalize();
    return 1;
  }

  char *mesh = argv[1];
  double tol = (argc > 2) ? atof(argv[2]) : 0.2;

  /* Read the geometry from the .re2 file */
  unsigned int nelt, nbcs;
  double *coord = NULL;
  long long *bcs = NULL;
  int nv;
  int err = parrsb_read_mesh(&nelt, &nv, NULL, &coord, &nbcs, &bcs, mesh,
                             MPI_COMM_WORLD, 1);
  check_error(err);

  /* Find connectivity */
  long long *vl = (long long *)calloc(nelt * nv, sizeof(long long));
  int ndim = nv == 8 ? 3 : 2;
  err |= parRSB_findConnectivity(vl, coord, nelt, ndim, bcs, nbcs, tol,
                                 MPI_COMM_WORLD, 0);
  check_error(err);

  /* Write connectivity file */
  err |= parrsb_dump_con(vl, nelt, nv, mesh, MPI_COMM_WORLD);
  check_error(err);

  if (vl != NULL)
    free(vl);
  if (coord != NULL)
    free(coord);
  if (bcs != NULL)
    free(bcs);

  MPI_Finalize();

  return 0;
}
