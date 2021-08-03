#include <math.h>
#include <stdio.h>

#include <mpi.h>

#include <gencon-impl.h>
#include <genmap.h>
#include <parRSB.h>

static int readRe2Header(Mesh *mesh_, MPI_File file, struct comm *c) {
  int rank = c->id;
  int size = c->np;
  MPI_Comm comm = c->c;

  char *buf = (char *)calloc(GC_RE2_HEADER_LEN + 1, sizeof(char));
  MPI_Status st;
  int err = MPI_File_read_all(file, buf, GC_RE2_HEADER_LEN, MPI_BYTE, &st);

  int nelgt, nelgv, nDim;
  char version[6];
  sscanf(buf, "%5s %d %d %d", version, &nelgt, &nDim, &nelgv);

  // TODO: Assert version
  int nVertex = (nDim == 2) ? 4 : 8;
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

  mesh_init(mesh_, nelt, nDim);
  Mesh mesh = *mesh_;
  mesh->nelgt = nelgt;
  mesh->nelgv = nelgv;
  mesh->nelt = nelt;
  mesh->nDim = nDim;

  free(buf);

  return 0;
}

static int readRe2Coordinates(Mesh mesh, MPI_File file, struct comm *c) {
  uint rank = c->id;
  uint size = c->np;
  MPI_Comm comm = c->c;

  int nelt = mesh->nelt;
  int nelgt = mesh->nelgt;
  int nDim = mesh->nDim;
  int nVertex = (nDim == 2) ? 4 : 8;

  slong out[2][1], bfr[2][1];
  slong in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0];

  int elemDataSize = nVertex * nDim * sizeof(double) + sizeof(double);
  int header_size = GC_RE2_HEADER_LEN + sizeof(float);

  /* calculate read size for element data on each MPI rank */
  int read_size = nelt * elemDataSize;
  if (rank == 0)
    read_size += header_size;

  char *buf = (char *)calloc(read_size, sizeof(char));
  char *buf0 = buf;
  MPI_Status st;
  int err = MPI_File_read_ordered(file, buf, read_size, MPI_BYTE, &st);
  if (err)
    return 1;

  if (rank == 0)
    buf0 += header_size;

  /* initialize array */
  uint nUnits = nelt * nVertex;
  array_init(struct Point_private, &mesh->elements, nUnits);
  Point ptr = mesh->elements.ptr;

  /* read elements for each rank */
  double x[GC_MAX_VERTICES], y[GC_MAX_VERTICES], z[GC_MAX_VERTICES];
  int i, j, k;
  if (nDim == 3) {
    for (i = 0; i < nelt; i++) {
      // skip group id
      buf0 += sizeof(double);
      READ_T(x, buf0, double, nVertex);
      buf0 += sizeof(double) * nVertex;
      READ_T(y, buf0, double, nVertex);
      buf0 += sizeof(double) * nVertex;
      READ_T(z, buf0, double, nVertex);
      buf0 += sizeof(double) * nVertex;

      for (k = 0; k < nVertex; k++) {
        ptr->x[0] = x[k];
        ptr->x[1] = y[k];
        ptr->x[2] = z[k];
        ptr->elementId = start + i;
        ptr->sequenceId = nVertex * (start + i) + k;
        ptr->origin = rank;
        ptr++;
      }
    }
  } else if (nDim == 2) {
    for (i = 0; i < nelt; i++) {
      // skip group id
      buf0 += sizeof(double);
      READ_T(x, buf0, double, nVertex);
      buf0 += sizeof(double) * nVertex;
      READ_T(y, buf0, double, nVertex);
      buf0 += sizeof(double) * nVertex;

      for (k = 0; k < nVertex; k++) {
        ptr->x[0] = x[k];
        ptr->x[1] = y[k];
        ptr->elementId = start + i;
        ptr->sequenceId = nVertex * (start + i) + k;
        ptr->origin = rank;
        ptr++;
      }
    }
  }

  mesh->elements.n = nUnits;

  free(buf);

  return 0;
}

static int readRe2Boundaries(Mesh mesh, MPI_File file, struct comm *c) {
  uint rank = c->id;
  uint size = c->np;
  MPI_Comm comm = c->c;

  int nelt = mesh->nelt;
  int nelgt = mesh->nelgt;
  int nDim = mesh->nDim;
  int nVertex = mesh->nVertex;

  int elemDataSize = nVertex * nDim * sizeof(double) + sizeof(double);
  int header_size = GC_RE2_HEADER_LEN + sizeof(float);

  MPI_Status st;
  char bufL[8];

  /* calculate offset for the curve side data */
  MPI_Offset curveOffset = header_size + nelgt * elemDataSize;
  if (rank == 0)
    MPI_File_read_at(file, curveOffset, bufL, sizeof(long), MPI_BYTE, &st);
  MPI_Bcast(bufL, sizeof(long), MPI_BYTE, 0, comm);

  double ncurvesD;
  READ_T(&ncurvesD, bufL, long, 1);
  long ncurves = ncurvesD;

  /* calculate offset for boundary conditions data */
  MPI_Offset boundaryOffset =
      curveOffset + sizeof(long) + sizeof(long) * 8 * ncurves;
  if (rank == 0)
    MPI_File_read_at(file, boundaryOffset, bufL, sizeof(long), MPI_BYTE, &st);
  MPI_Bcast(bufL, sizeof(long), MPI_BYTE, 0, comm);

  double nbcsD;
  READ_T(&nbcsD, bufL, long, 1);
  long nbcs = nbcsD;

  int nbcsLocal = nbcs / size, nrem = nbcs - nbcsLocal * size;
  nbcsLocal += (rank >= (size - nrem) ? 1 : 0);

  slong out[2][1], bfr[2][1], in = nbcsLocal;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0];

  int offset = boundaryOffset + sizeof(long) + start * 8 * sizeof(long);
  int read_size = nbcsLocal * sizeof(long) * 8;
  char *buf = calloc(read_size, sizeof(char)), *buf0 = buf;
  MPI_File_read_at_all(file, offset, buf, read_size, MPI_BYTE, &st);

  double tmp[5];
  char cbc[4];
  struct Boundary_private boundary;
  sint i;
  for (i = 0; i < nbcsLocal; i++) {
    READ_T(tmp, buf0, long, 1);
    buf0 += sizeof(long);
    boundary.elementId = (long)tmp[0];

    READ_T(tmp, buf0, long, 1);
    buf0 += sizeof(long);
    boundary.faceId = (long)tmp[0];

    READ_T(tmp, buf0, long, 5);
    buf0 += 5 * sizeof(long);
    READ_T(cbc, buf0, char, 3);
    buf0 += sizeof(long);
    cbc[3] = '\0';

    if (strcmp(cbc, GC_PERIODIC) == 0) {
      boundary.bc[0] = (long)tmp[0];
      boundary.bc[1] = (long)tmp[1];
      array_cat(struct Boundary_private, &mesh->boundary, &boundary, 1);
    }
  }
  free(buf);
}

static int read_geometry(Mesh *mesh_, char *fname, struct comm *c) {
  int nelt, nDim, nVertex;
  int errs = 0;

  uint rank = c->id;
  uint size = c->np;
  MPI_Comm comm = c->c;

  MPI_File file;
  int err = MPI_File_open(comm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
  if (err != 0) {
    if (rank == 0)
      printf("%s:%d Error opening file: %s\n", __FILE__, __LINE__, fname);
    return 1;
  }

  readRe2Header(mesh_, file, c);
  Mesh mesh = *mesh_;
  readRe2Coordinates(mesh, file, c);
  readRe2Boundaries(mesh, file, c);

  err = MPI_File_close(&file);
  if (err)
    errs++;

  MPI_Barrier(comm);

  return errs;
}

static int read_connectivity(Mesh mesh, char *fname, struct comm *c) {
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

int parrsb_read_mesh(unsigned int *nel, int *nv, long long **vl, double **coord,
                     unsigned int *nbcs, long long **bcs, char *name,
                     MPI_Comm comm, int read) {
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
    get_bcs(nbcs, bcs, mesh);
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
  MPI_Bcast(&ndim, 1, MPI_INT, 0, MPI_COMM_WORLD);
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

  int i, j, temp;
  for (i = 0; i < nelt; i++) {
    temp = start + i + 1;
    WRITE_INT(buf0, temp);
    buf0 += sizeof(int);
    for (j = 0; j < nv; j++) {
      temp = vl[i * nv + j];
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
