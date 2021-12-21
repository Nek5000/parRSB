#include <stdio.h>

#include <genmap-impl.h>

#define write_T(dest, val, T, nunits)                                          \
  do {                                                                         \
    memcpy(dest, (val), sizeof(T) * nunits);                                   \
    dest += sizeof(T) * nunits;                                                \
  } while (0)

#define write_1(dest, val, T)                                                  \
  do {                                                                         \
    memcpy(dest, &(val), sizeof(T));                                           \
  } while (0)

int GenmapFiedlerDump(const char *fname, struct rsb_element *elm, uint nelt,
                      int nv, struct comm *c) {
  MPI_File file;
  int err = MPI_File_open(c->c, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &file);

  uint rank = c->id;
  if (err != 0 && rank == 0)
    fprintf(stderr, "%s:%d Error opening file %s for writing.\n", __FILE__,
            __LINE__, fname);
  if (err != 0)
    return err;

  slong out[2][1], buf[2][1];
  slong in = nelt;
  comm_scan(out, c, gs_long, gs_add, &nelt, 1, buf);
  slong start = out[0][0];
  slong nelgt = out[1][0];

  int ndim = (nv == 8) ? 3 : 2;
  uint write_size =
      ((ndim + 1) * sizeof(double) + sizeof(GenmapLong) + sizeof(uint)) * nelt;
  if (rank == 0)
    write_size += sizeof(long) + sizeof(int); // for nelgt and ndim

  char *pbuf, *pbuf0;
  pbuf = pbuf0 = (char *)calloc(write_size, sizeof(char));
  if (rank == 0) {
    write_T(pbuf0, &nelgt, slong, 1);
    write_T(pbuf0, &ndim, int, 1);
  }

  uint i;
  for (i = 0; i < nelt; i++) {
    write_T(pbuf0, &elm[i].globalId, GenmapULong, 1);
    write_T(pbuf0, elm[i].coord, double, ndim);
    write_T(pbuf0, &elm[i].fiedler, double, 1);
    write_T(pbuf0, &rank, uint, 1);
  }

  MPI_Status st;
  err = MPI_File_write_ordered(file, pbuf, write_size, MPI_BYTE, &st);
  if (err != 0 && rank == 0)
    fprintf(stderr, "%s:%d Error opening file %s for writing.\n", __FILE__,
            __LINE__, fname);

  err += MPI_File_close(&file);
  genmap_barrier(c);

  free(pbuf);

  return err;
}

int GenmapElementDump(const char *fname, struct rsb_element *elm, uint nelt,
                      int nv, struct comm *c, int dump) {
  if (dump == 0)
    return 1;

  MPI_File file;
  int err = MPI_File_open(c->c, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &file);

  uint rank = c->id;
  if (err != 0 && rank == 0)
    fprintf(stderr, "%s:%d Error opening file %s for writing.\n", __FILE__,
            __LINE__, fname);
  if (err != 0)
    return err;

  slong out[2][1], buf[2][1];
  slong in = nelt;
  comm_scan(out, c, gs_long, gs_add, &nelt, 1, buf);
  slong start = out[0][0];
  slong nelgt = out[1][0];

  int ndim = (nv == 8) ? 3 : 2;
  uint write_size = sizeof(int) * nelt;
  if (rank == 0)
    write_size += sizeof(int) + sizeof(int); // for nelgt and ndim

  if (c->id == 0)
    printf("nelgt = %d slong=%u int=%u ulong=%u\n", nelgt, sizeof(slong),
           sizeof(int), sizeof(GenmapULong));

  char *pbuf, *pbuf0;
  pbuf = pbuf0 = (char *)calloc(write_size, sizeof(char));
  int temp;
  if (rank == 0) {
    temp = nelgt;
    write_1(pbuf0, temp, int);
    pbuf0 += sizeof(int);
    write_1(pbuf0, ndim, int);
    pbuf0 += sizeof(int);
  }

  uint i;
  for (i = 0; i < nelt; i++) {
    temp = elm[i].globalId;
    write_1(pbuf0, temp, int);
    pbuf0 += sizeof(int);
  }

  MPI_Status st;
  err = MPI_File_write_ordered(file, pbuf, write_size, MPI_BYTE, &st);
  if (err != 0 && rank == 0)
    fprintf(stderr, "%s:%d Error opening file %s for writing.\n", __FILE__,
            __LINE__, fname);

  err += MPI_File_close(&file);
  genmap_barrier(c);

  free(pbuf);

  return err;
}

int GenmapVectorDump(const char *fname, GenmapScalar *y,
                     struct rsb_element *elm, uint nelt, int nv,
                     struct comm *c) {
  MPI_File file;
  int err = MPI_File_open(c->c, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &file);

  uint rank = c->id;
  if (err != 0 && rank == 0)
    fprintf(stderr, "%s:%d Error opening file %s for writing.\n", __FILE__,
            __LINE__, fname);
  if (err != 0)
    return err;

  slong out[2][1], buf[2][1];
  slong in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
  slong start = out[0][0];
  slong nelgt = out[1][0];

  int ndim = (nv == 8) ? 3 : 2;
  uint write_size = ((ndim + 1) * sizeof(double) + sizeof(GenmapLong)) * nelt;
  if (rank == 0)
    write_size += sizeof(long) + sizeof(int); // for nelgt and ndim

  char *pbuf, *pbuf0;
  pbuf = pbuf0 = (char *)calloc(write_size, sizeof(char));
  if (rank == 0) {
    write_T(pbuf0, &nelgt, slong, 1);
    write_T(pbuf0, &ndim, int, 1);
  }

  uint i;
  for (i = 0; i < nelt; i++) {
    write_T(pbuf0, &elm[i].globalId, GenmapULong, 1);
    write_T(pbuf0, &elm[i].coord[0], double, ndim);
    write_T(pbuf0, &y[i], double, 1);
  }

  MPI_Status st;
  err = MPI_File_write_ordered(file, pbuf, write_size, MPI_BYTE, &st);
  if (err != 0 && rank == 0)
    fprintf(stderr, "%s:%d Error opening file %s for writing.\n", __FILE__,
            __LINE__, fname);

  err += MPI_File_close(&file);
  genmap_barrier(c);

  free(pbuf);

  return err;
}

#undef write_T
#undef write_1
