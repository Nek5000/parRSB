//=============================================================================
// Test ILU factorization
//
#include "ilu.h"
#include "parRSB.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  struct parrsb_input *in = parrsb_parse_input(argc, argv, comm);
  int err = (in == NULL);
  parrsb_check_error(err, comm);

  int id;
  MPI_Comm_rank(comm, &id);

  // Read the geometry from the .re2 file, find connectiviy, partition and then
  // distribute the mesh.
  unsigned int nelt;
  int nv;
  long long *vl = NULL;
  double *coord = NULL;
  parrsb_setup_mesh(&nelt, &nv, &vl, &coord, in, comm);

  // Setup ILU
  ilu_options iluopt = {.type = in->type,
                        .verbose = in->verbose,
                        .tol = in->tol,
                        .pivot = 0,
                        .nnz_per_row = 0};
  struct ilu *ilu = ilu_setup(nelt, nv, vl, &iluopt, comm);
  ilu_free(ilu);

  // Free resources
  free(vl), free(coord);
  free(in);
  MPI_Finalize();

  return 0;
}
