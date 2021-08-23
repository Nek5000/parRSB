/*
 * Generate partitions (.ma2) from Nek5000's mesh (.re2) file.
 */
#include <parRSB.h>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  parrsb_input *in = parrsb_parse_input(argc, argv);

  int id;
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  int color = 0;
  if (id < in->nactive)
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
    err = parrsb_read_mesh(&nelt, &nv, NULL, &coord, &nbcs, &bcs, in->mesh,
                           comm, 1);
  parrsb_check_error(err, comm);

  /* Find connectivity */
  long long *vl = (long long *)calloc(nelt * nv, sizeof(long long));
  int ndim = nv == 8 ? 3 : 2;
  if (color == 1)
    err |= parRSB_findConnectivity(vl, coord, nelt, ndim, bcs, nbcs, in->tol,
                                   comm, 0);
  parrsb_check_error(err, comm);

  if (color == 1)
    parrsb_print_part_stat(vl, nelt, nv, comm);

  /* Partition the mesh */
  parrsb_options options = parrsb_default_options;
  int *part = (int *)calloc(nelt, sizeof(int));
  if (color == 1)
    err |= parRSB_partMesh(part, NULL, vl, coord, nelt, nv, options, comm);
  parrsb_check_error(err, comm);

  /* Redistribute data */
  if (color == 1)
    err |= parrsb_distribute_elements(&nelt, &vl, &coord, part, nv, comm);
  parrsb_check_error(err, comm);

  if (color == 1)
    parrsb_print_part_stat(vl, nelt, nv, comm);

  /* Write map file */
  if (color == 1 && in->dump == 1)
    err |= parrsb_dump_map(nelt, nv, part, vl, in->mesh, comm);
  parrsb_check_error(err, comm);

  /* Free resources */
  if (part != NULL)
    free(part);
  if (vl != NULL)
    free(vl);
  if (coord != NULL)
    free(coord);
  if (bcs != NULL)
    free(bcs);
  if (in != NULL)
    free(in);

  MPI_Comm_free(&comm);
  MPI_Finalize();

  return 0;
}
