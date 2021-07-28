/*
 * Generate partitions (.ma2) from Nek5000's mesh (.re2) file.
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
  MPI_Comm_split(MPI_COMM_WORLD, color, id, &comm);

  /* Read the geometry from the .re2 file */
  unsigned int nelt, nbcs;
  double *coord = NULL;
  long long *bcs = NULL;
  int nv;
  int err = 0;
  if (color == 1)
    err =
        parrsb_read_mesh(&nelt, &nv, NULL, &coord, &nbcs, &bcs, mesh, comm, 1);
  check_error(err);

  /* Find connectivity */
  long long *vl = (long long *)calloc(nelt * nv, sizeof(long long));
  int ndim = nv == 8 ? 3 : 2;
  if (color == 1)
    err |=
        parRSB_findConnectivity(vl, coord, nelt, ndim, bcs, nbcs, tol, comm, 0);
  check_error(err);

  if (color == 1)
    parrsb_part_stat(vl, nelt, nv, comm);

  /* Partition the mesh */
  parRSB_options options = parrsb_default_options;
  int *part = (int *)calloc(nelt, sizeof(int));
  if (color == 1)
    err |= parRSB_partMesh(part, NULL, vl, coord, nelt, nv, &options, comm);
  check_error(err);

  /* Redistribute data */
  if (color == 1)
    err |= parrsb_distribute_elements(&nelt, &vl, &coord, part, nv, comm);
  check_error(err);

  if (color == 1)
    parrsb_part_stat(vl, nelt, nv, comm);

  /* Write map file */
  if (color == 1)
    err |= parrsb_dump_map(nelt, nv, part, vl, mesh, comm);
  check_error(err);

  /* Free resources */
  if (part != NULL)
    free(part);
  if (vl != NULL)
    free(vl);
  if (coord != NULL)
    free(coord);
  if (bcs != NULL)
    free(bcs);

  MPI_Comm_free(&comm);
  MPI_Finalize();

  return 0;
}
