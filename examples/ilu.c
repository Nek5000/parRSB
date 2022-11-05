//=============================================================================
// Test ILU factorization
//
#include "ilu.h"
#include "parRSB.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  struct parrsb_cmd_opts *in = parrsb_parse_cmd_opts(argc, argv);
  int err = (in == NULL);
  parrsb_check_error(err, comm);

  // Read the geometry from the .re2 file, find connectiviy, partition and then
  // distribute the mesh.
  unsigned int ne, nv;
  long long *vl = NULL;
  double *coord = NULL;
  parrsb_setup_mesh(&ne, &nv, &vl, &coord, in, comm);

  buffer bfr;
  buffer_init(&bfr, 1024);

  // Setup ILU
  ilu_options opts = {.type = in->ilu_type,
                      .tol = in->ilu_tol,
                      .pivot = in->ilu_pivot,
                      .verbose = in->verbose,
                      .nnz_per_row = in->ilu_nnz_per_row};
  struct ilu *ilu = ilu_setup(ne, nv, vl, &opts, comm, &bfr);

  if (in->ilu_solve) {
    double *x = NULL, *b = NULL;
    ilu_solve(x, ilu, b, &bfr);
  }

  ilu_free(ilu), buffer_free(&bfr);
  free(vl), free(coord), free(in);
  MPI_Finalize();

  return 0;
}
