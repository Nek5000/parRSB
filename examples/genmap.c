/*
 * Parition mesh using Nek5000's connectivity (.co2) and map (.ma2) file.
 */
#include <mpi.h>
#include <stdio.h>

#include <parRSB.h>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int id, np;
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  if (argc < 2 || argc > 4) {
    if (id == 0)
      printf("Usage: %s <mesh file> [ranks] [tol]\n", argv[0]);
    MPI_Finalize();
    return 1;
  }

  char *mesh = argv[1];
  double tol = (argc > 2) ? atof(argv[2]) : 0.2;
  int ranks = (argc > 3) ? atoi(argv[3]) : np;

  int color = 0;
  if (id < ranks)
    color = 1;
  MPI_Comm comm;
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &comm);

  /* Read the geometry from the .re2 file */
  unsigned int nelt;
  int nv;
  double *coord = NULL;
  int err = parrsb_read_mesh(&nelt, &nv, NULL, &coord, mesh, comm, 1);

  /* Find connectivity */
  long long *vl = NULL;
  if (err == 0) {
    vl = (long long *)calloc(nelt * nv, sizeof(long long));
    int ndim = nv == 8 ? 3 : 2;
    err |=
        parRSB_findConnectivity(vl, coord, nelt, ndim, NULL, 0, tol, comm, 0);
  }

  parrsb_part_stat(vl, nelt, nv, comm);

  /* Partition the mesh */
  parRSB_options options = parrsb_default_options;
  int *part = (int *)calloc(nelt, sizeof(int));
  err |= parRSB_partMesh(part, NULL, vl, coord, nelt, nv, &options, comm);

  /* Redistribute data */
  err |= parrsb_distribute_elements(&nelt, &vl, &coord, part, nv, comm);

  parrsb_part_stat(vl, nelt, nv, MPI_COMM_WORLD);

  /* Write map file */

  if (part != NULL)
    free(part);
  if (vl != NULL)
    free(vl);
  if (coord != NULL)
    free(coord);

  MPI_Comm_free(&comm);

  MPI_Finalize();

  return 0;
}
