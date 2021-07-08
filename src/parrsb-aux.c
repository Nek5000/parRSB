#include <mpi.h>
#include <stdio.h>

#include <gencon.h>
#include <genmap.h>

#define MAXNV 8  /* maximum number of vertices per element */
#define MAXDIM 3 /* maximum number of vertices per element */

/* elem_data must always start with vtx_data */
typedef struct {
  int proc;
  long long vtx[MAXNV];
} vtx_data;

typedef struct {
  int proc;
  long long vtx[MAXNV];
  double coord[MAXNV * MAXDIM];
} elem_data;

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

int parrsb_distribute_elements(unsigned int *nelt_, long long **vl_,
                               double **coord_, int *part, int nv,
                               MPI_Comm comm) {
  int ndim = (nv == 8) ? 3 : 2;
  long long *vl = *vl_;
  double *coord = *coord_;

  size_t unit_size = 0;
  if (coord != NULL)
    unit_size = sizeof(elem_data);
  else
    unit_size = sizeof(vtx_data);

  uint nelt = *nelt_;
  struct array elements;
  array_init_(&elements, nelt, unit_size, __FILE__, __LINE__);

  elem_data data;
  int e, n;
  for (e = 0; e < nelt; ++e) {
    data.proc = part[e];
    for (n = 0; n < nv; ++n)
      data.vtx[n] = vl[e * nv + n];
    array_cat_(unit_size, &elements, &data, 1, __FILE__, __LINE__);
  }
  assert(elements.n == nelt);

  if (coord != NULL) {
    elem_data *ed = elements.ptr;
    for (e = 0; e < nelt; e++)
      for (n = 0; n < ndim * nv; n++)
        ed[e].coord[n] = coord[e * ndim * nv + n];
  }

  struct comm c;
  comm_init(&c, comm);

  struct crystal cr;
  crystal_init(&cr, &c);

  sarray_transfer_(&elements, unit_size, offsetof(vtx_data, proc), 0, &cr);

  nelt = elements.n;
  *nelt_ = nelt;

  vl = *vl_ = (long long *)realloc(*vl_, nv * nelt * sizeof(long long));
  for (e = 0; e < nelt; ++e) {
    vtx_data *vd = (vtx_data *)(elements.ptr + unit_size * e);
    for (n = 0; n < nv; ++n)
      vl[e * nv + n] = vd->vtx[n];
  }

  if (coord != NULL) {
    coord = *coord_ =
        (double *)realloc(*coord_, ndim * nv * nelt * sizeof(double));
    elem_data *ed = elements.ptr;
    for (e = 0; e < nelt; ++e) {
      for (n = 0; n < ndim * nv; ++n)
        coord[e * ndim * nv + n] = ed[e].coord[n];
    }
  }

  crystal_free(&cr);
  comm_free(&c);
  array_free(&elements);

  return 0;
}

#define WRITE_INT(dest, val)                                                   \
  do {                                                                         \
    memcpy(dest, &(val), sizeof(int));                                         \
  } while (0)

int parrsb_dump_con(long long *vl, unsigned int nelt, int nv, char *name,
                    MPI_Comm comm) {
  const char version[5] = "#v001";
  const float test = 6.54321;

  struct comm c;
  comm_init(&c, comm);

  uint id = c.id;

  MPI_File file;
  int err = MPI_File_open(comm, name, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &file);
  if (err) {
    if (id == 0) {
      fprintf(stderr, "%s:%d Error opening file: %s for writing.\n", __FILE__,
              __LINE__, name);
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
#if 0
    printf("%5s%12d%12d%12d\n", version, (int)nelgt, (int)nelgv, nv);
#endif
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

void parrsb_part_stat(long long *vtx, int nel, int nv, MPI_Comm ce) {
  int i, j;

  struct comm comm;
  int np, id;

  int Nmsg;
  int *Ncomm;

  int nelMin, nelMax;
  int ncMin, ncMax, ncSum;
  int nsMin, nsMax, nsSum;
  int nssMin, nssMax, nssSum;

  struct gs_data *gsh;
  int b;

  int numPoints;
  long long *data;

  comm_init(&comm, ce);
  np = comm.np;
  id = comm.id;

  if (np == 1)
    return;

  numPoints = nel * nv;
  data = (long long *)malloc((numPoints + 1) * sizeof(long long));
  for (i = 0; i < numPoints; i++)
    data[i] = vtx[i];

  gsh = gs_setup(data, numPoints, &comm, 0, gs_pairwise, 0);

  pw_data_nmsg(gsh, &Nmsg);
  Ncomm = (int *)malloc((Nmsg + 1) * sizeof(int));
  pw_data_size(gsh, Ncomm);

  gs_free(gsh);
  free(data);

  ncMax = Nmsg;
  ncMin = Nmsg;
  ncSum = Nmsg;
  comm_allreduce(&comm, gs_int, gs_max, &ncMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &ncMin, 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &ncSum, 1, &b);

  nsMax = Ncomm[0];
  nsMin = Ncomm[0];
  nsSum = Ncomm[0];
  for (i = 1; i < Nmsg; ++i) {
    nsMax = Ncomm[i] > Ncomm[i - 1] ? Ncomm[i] : Ncomm[i - 1];
    nsMin = Ncomm[i] < Ncomm[i - 1] ? Ncomm[i] : Ncomm[i - 1];
    nsSum += Ncomm[i];
  }
  comm_allreduce(&comm, gs_int, gs_max, &nsMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nsMin, 1, &b);

  nssMin = nsSum;
  nssMax = nsSum;
  nssSum = nsSum;
  comm_allreduce(&comm, gs_int, gs_max, &nssMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nssMin, 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &nssSum, 1, &b);

  if (Nmsg)
    nsSum = nsSum / Nmsg;
  else
    nsSum = 0;
  comm_allreduce(&comm, gs_int, gs_add, &nsSum, 1, &b);

  nelMax = nel;
  nelMin = nel;
  comm_allreduce(&comm, gs_int, gs_max, &nelMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nelMin, 1, &b);

  if (id == 0) {
    printf(" Max neighbors: %d | Min neighbors: %d | Avg neighbors: %lf\n",
           ncMax, ncMin, (double)ncSum / np);
    printf(" Max nvolume: %d | Min nvolume: %d | Avg nvolume: %lf\n", nsMax,
           nsMin, (double)nsSum / np);
    printf(" Max volume: %d | Min volume: %d | Avg volume: %lf\n", nssMax,
           nssMin, (double)nssSum / np);
    printf(" Max elements: %d | Min elements: %d\n", nelMax, nelMin);
    fflush(stdout);
  }

  free(Ncomm);
  comm_free(&comm);
}

#undef MAXNV
