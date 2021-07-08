#include <mpi.h>
#include <stdio.h>

#include <gencon.h>
#include <genmap.h>
#include <parRSB.h>

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
