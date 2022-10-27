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
  ilu_options iluopt = {.type = in->ilu_type,
                        .tol = in->ilu_tol,
                        .null_space = in->ilu_null_space,
                        .pivot = in->ilu_pivot,
                        .verbose = in->verbose,
                        .nnz_per_row = 0};
  struct ilu *ilu = ilu_setup(ne, nv, vl, &iluopt, comm, &bfr);

  struct comm c;
  comm_init(&c, comm);
  slong out[2][1], wrk[2][1], n = ne;
  comm_scan(out, &c, gs_long, gs_add, &n, 1, wrk);
  ulong s = out[0][0];
  comm_free(&c);

  double *x = tcalloc(double, 2 * ne), *b = x + ne;
  for (unsigned i = 0; i < ne; i++)
    b[i] = s + i + 1;

  // ilu_solve(x, ilu, b, &bfr);

  ilu_free(ilu), buffer_free(&bfr);
  free(x), free(vl), free(coord), free(in);
  MPI_Finalize();

  return 0;
}
