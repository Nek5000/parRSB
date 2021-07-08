#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include <parRSB.h>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int id;
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  if (argc < 2) {
    if (id == 0)
      printf("Usage: ./%s <mesh file> [tol]\n", argv[0]);
    MPI_Finalize();
    return 1;
  }

  char *mesh = argv[1];
  double tol = (argc > 2) ? atof(argv[2]) : 0.2;

  /* Read the geometry from the .re2 file */
  unsigned int nelt;
  int nv;
  double *coord = NULL;
  int err = parrsb_read_mesh(&nelt, &nv, NULL, &coord, mesh, MPI_COMM_WORLD, 1);

  long long *vl = NULL;
  if (err == 0) {
    vl = (long long *)calloc(nelt * nv, sizeof(long long));
    int ndim = nv == 8 ? 3 : 2;
    err |= parRSB_findConnectivity(vl, coord, nelt, ndim, NULL, 0, tol,
                                   MPI_COMM_WORLD, 0);
  }

  /* Write connectivity */

  if (vl != NULL)
    free(vl);
  if (coord != NULL)
    free(coord);

  MPI_Finalize();

  return 0;
}
