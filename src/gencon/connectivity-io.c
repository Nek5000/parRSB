#include <gencon-impl.h>
#include <genmap-impl.h>

int read_connectivity(Mesh mesh, char *fname, struct comm *c) {
  comm_ext comm = c->c;
  int rank = c->id;
  int size = c->np;

  MPI_File file;
  int err = MPI_File_open(comm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
  if (err != 0) {
    if (rank == 0)
      printf("Error opening %s for reading.\n", fname);
    return 1;
  }

  char *buf = (char *)calloc(GC_CO2_HEADER_LEN + 1, sizeof(char));
  MPI_Status st;
  err = MPI_File_read_all(file, buf, GC_CO2_HEADER_LEN, MPI_BYTE, &st);

  int nelgt, nelgv, nVertex;
  char version[6];
  sscanf(buf, "%5s%12d%12d%12d", version, &nelgt, &nelgv, &nVertex);

  /* TODO: Assert version */
  int nDim = (nVertex == 8) ? 3 : 2;
  int nelt = nelgt / size;
  int nrem = nelgt - nelt * size;
  nelt += (rank > (size - 1 - nrem) ? 1 : 0);

  float byte_test;
  MPI_File_read_all(file, &byte_test, 4, MPI_BYTE, &st);
  if (fabs(byte_test - 6.543210) > 1e-7) {
    if (rank == 0)
      printf("ERROR byte_test failed! %f\n", byte_test);
    return 1;
  }

  // Check for the same .re2 mesh
  assert(mesh->nDim == nDim);
  assert(mesh->nelgt == nelgt);
  assert(mesh->nelgv == nelgv);
  assert(mesh->nelt == nelt);
  assert(mesh->elements.n == nelt * nVertex);

#if 0
  printf("co2 %d: %s %d %d %d %d\n", rank, version, nelgt, nelgv, nelt, nDim);
#endif

  slong out[2][1], bfr[2][1], in[1];
  in[0] = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0];

  int read_size = nelt * (nVertex + 1) * sizeof(int);
  int header_size = GC_CO2_HEADER_LEN + sizeof(float);
  if (rank == 0)
    read_size += header_size;

  buf = (char *)realloc(buf, read_size * sizeof(char));
  err = MPI_File_read_ordered(file, buf, read_size, MPI_BYTE, &st);
  err = MPI_File_close(&file);

  char *buf0 = buf;
  if (rank == 0)
    buf0 += header_size;

  Point ptr = mesh->elements.ptr;
  int i, j, tmp1, tmp2;
  for (i = 0; i < nelt; i++) {
    READ_T(&tmp1, buf0, int, 1);
    buf0 += sizeof(int);
    for (j = 0; j < nVertex; j++) {
      ptr[i * nVertex + j].elementId = tmp1;
      READ_T(&tmp2, buf0, int, 1);
      buf0 += sizeof(int);
      ptr[i * nVertex + j].globalId = tmp2;
    }
  }

  free(buf);

  return 0;
}
