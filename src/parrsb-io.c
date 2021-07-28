#include <math.h>
#include <stdio.h>

#include <mpi.h>

#include <gencon.h>
#include <genmap.h>
#include <parRSB.h>

int parrsb_read_mesh(unsigned int *nel, int *nv, long long **vl, double **coord,
                     char *name, MPI_Comm comm, int read) {
  struct comm c;
  comm_init(&c, comm);

  Mesh mesh;

  /* Read geometry from .re2 file */
  if (read & 1) {
    char geom_name[BUFSIZ];
    strncpy(geom_name, name, BUFSIZ);
    strncat(geom_name, ".re2", 5);

    read_geometry(&mesh, geom_name, &c);
    get_vertex_coordinates(coord, mesh);
  }

  /* Read connectivity from .co2 file */
  if (read & 2) {
    char conn_name[BUFSIZ];
    strncpy(conn_name, name, BUFSIZ);
    strncat(conn_name, ".co2", 5);

    assert(read & 1);
    read_connectivity(mesh, conn_name, &c);
    get_vertex_ids(vl, mesh);
  }

  *nel = (unsigned int)get_mesh_nel(mesh);

  int ndim = get_mesh_dim(mesh);
#if defined(MPI)
  MPI_Bcast(&ndim, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  *nv = (ndim == 3) ? 8 : 4;

  mesh_free(mesh);
  comm_free(&c);

  return 0;
}

#define WRITE_INT(dest, val)                                                   \
  do {                                                                         \
    memcpy(dest, &(val), sizeof(int));                                         \
  } while (0)

int parrsb_dump_con(long long *vl, unsigned int nelt, int nv, char *name,
                    MPI_Comm comm) {
  const char version[6] = "#v001";
  const float test = 6.54321;

  struct comm c;
  comm_init(&c, comm);
  uint id = c.id;

  char co2_name[BUFSIZ];
  strncpy(co2_name, name, BUFSIZ);
  strncat(co2_name, ".co2", 5);

  MPI_File file;
  int err = MPI_File_open(comm, co2_name, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &file);
  if (err) {
    if (id == 0) {
      fprintf(stderr, "%s:%d Error opening file: %s for writing.\n", __FILE__,
              __LINE__, co2_name);
      fflush(stdout);
    }
    return err;
  }

  slong out[2][1], bfr[2][1];
  slong in = nelt;
  comm_scan(out, &c, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0];
  slong nelgt = out[1][0];
  slong nelgv = nelgt;

  int write_size = nelt * (nv + 1) * sizeof(int);
  int header_size = GC_CO2_HEADER_LEN + sizeof(float);
  if (id == 0)
    write_size += header_size;

  char *buf = (char *)calloc(write_size, sizeof(char));
  char *buf0 = buf;

  if (id == 0) {
    sprintf(buf0, "%5s%12d%12d%12d", version, (int)nelgt, (int)nelgv, nv);
    memset(buf0 + strlen(buf0), ' ', GC_CO2_HEADER_LEN - strlen(buf0));
    buf0[GC_CO2_HEADER_LEN] = '\0';
    buf0 += GC_CO2_HEADER_LEN;
    memcpy(buf0, &test, sizeof(float));
    buf0 += sizeof(float);
  }

  int i, j, k, temp;
  for (i = 0; i < nelt; i++) {
    temp = start + i + 1;
    WRITE_INT(buf0, temp);
    buf0 += sizeof(int);
    for (j = 0; j < nv; j++) {
      k = PRE_TO_SYM_VERTEX[j];
      temp = vl[i * nv + k];
      WRITE_INT(buf0, temp);
      buf0 += sizeof(int);
    }
  }

  MPI_Status st;
  err |= MPI_File_write_ordered(file, buf, write_size, MPI_BYTE, &st);
  err |= MPI_File_close(&file);

  comm_barrier(&c);
  comm_free(&c);

  free(buf);

  return err;
}

#undef WRITE_INT

#define HEADER_LEN 132

int parrsb_dump_map(int nelt, int nv, int *pmap, long long *vtx, char *fileName,
                    MPI_Comm comm) {
  char version[6] = "#v001";
  float test = 6.54321;

  int nelgt = nelt;
  MPI_Allreduce(&nelt, &nelgt, 1, MPI_INT, MPI_SUM, comm);

  const int npts = nelgt * nv;
  const int depth = (int)log2(1.0 * nelgt);
  const int d2 = (int)(pow(2, depth) + 0.5);
  int nactive = nelgt;
  int nrnk = nelgt;
  int noutflow = 0;

  char header[HEADER_LEN];
  header[HEADER_LEN] = '\0';
  sprintf(header, "%5s%12d%12d%12d%12d%12d%12d%12d", version, nelgt, nactive,
          depth, d2, npts, nrnk, noutflow);
  memset(header + strlen(header), ' ', HEADER_LEN - strlen(header));
  header[HEADER_LEN] = '\0';

  MPI_Info infoIn;
  MPI_Info_create(&infoIn);
  MPI_Info_set(infoIn, "access_style", "write_once,random");

  int errs = 0;
  MPI_File file;
  int err = MPI_File_open(comm, fileName, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                          infoIn, &file);
  if (err) {
    errs++;
    MPI_Abort(comm, 911);
  }

  int rank;
  MPI_Comm_rank(comm, &rank);

  int writeSize = 0;
  if (rank == 0)
    writeSize = HEADER_LEN * sizeof(char) + sizeof(float);
  writeSize += (nv + 1) * nelt * sizeof(int);

  char *buf = (char *)malloc(writeSize), *buf0 = buf;

  if (rank == 0) {
    memcpy(buf0, header, HEADER_LEN * sizeof(char)), buf0 += HEADER_LEN;
    memcpy(buf0, &test, sizeof(float)), buf0 += sizeof(float);
  }

  int ivtx[8];
  int i, j;
  for (i = 0; i < nelt; i++) {
    memcpy(buf0, &pmap[i], sizeof(int));
    buf0 += sizeof(int);

    for (j = 0; j < nv; j++)
      ivtx[j] = vtx[i * nv + j];
    memcpy(buf0, ivtx, sizeof(int) * nv);
    buf0 += nv * sizeof(int);
  }

  MPI_Status status;
  err = MPI_File_write_ordered(file, buf, writeSize, MPI_BYTE, &status);
  if (err)
    errs++;

  err = MPI_File_close(&file);
  if (err)
    errs++;

  MPI_Barrier(comm);
  free(buf);

  return errs;
}

#undef HEADER_LEN
