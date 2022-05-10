//=============================================================================
// Test ILU factorization
//
#include <parRSB.h>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  struct parrsb_input *in = parrsb_parse_input(argc, argv);
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
  if (vl == NULL)
    err = 1;
  parrsb_check_error(err, MPI_COMM_WORLD);

  int ndim = nv == 8 ? 3 : 2;
  err = parrsb_find_conn(vl, coord, nelt, ndim, bcs, nbcs, in->tol,
                         MPI_COMM_WORLD, 0);
  parrsb_check_error(err, MPI_COMM_WORLD);

  // Partition the mesh
  int *part = (int *)calloc(nelt, sizeof(int));
  if (part == NULL)
    err = 1;
  parrsb_check_error(err, MPI_COMM_WORLD);

  parrsb_options options = parrsb_default_options;
  err = parrsb_part_mesh(part, NULL, vl, coord, nelt, nv, options,
                         MPI_COMM_WORLD);
  parrsb_check_error(err, MPI_COMM_WORLD);

  // Redistribute data based on identified partitions
  err =
      parrsb_distribute_elements(&nelt, &vl, &coord, part, nv, MPI_COMM_WORLD);
  parrsb_check_error(err, MPI_COMM_WORLD);

  // Setup ILU
  struct comm c;
  comm_init(&c, MPI_COMM_WORLD);

#if 0
  struct ilu *ilu = ilu_setup(nelt, nv, vl, 0, 0, 0, &c, 1);
  ilu_free(ilu);
#endif

  comm_free(&c);

  return 0;
}
