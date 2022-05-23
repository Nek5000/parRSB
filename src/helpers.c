#include <getopt.h>
#include <stdio.h>

#include <genmap-impl.h>
#include <parRSB.h>

int parrsb_distribute_elements(unsigned int *nelt_, long long **vl_,
                               double **coord_, int *part, int nv,
                               MPI_Comm comm) {
  typedef struct {
    int proc;
    long long vtx[MAXNV];
    double coord[MAXNV * MAXDIM];
  } elem_data;

  uint nelt = *nelt_;
  struct array elements;
  array_init(elem_data, &elements, nelt);

  elem_data data;
  long long *vl = *vl_;
  uint e, n;
  for (e = 0; e < nelt; ++e) {
    data.proc = part[e];
    for (n = 0; n < nv; ++n)
      data.vtx[n] = vl[e * nv + n];
    array_cat(elem_data, &elements, &data, 1);
  }
  assert(elements.n == nelt);

  int ndim = (nv == 8) ? 3 : 2;
  elem_data *ed = elements.ptr;
  double *coord = (coord_ == NULL ? NULL : *coord_);
  if (coord != NULL) {
    for (e = 0; e < nelt; e++)
      for (n = 0; n < ndim * nv; n++)
        ed[e].coord[n] = coord[e * ndim * nv + n];
  }

  struct comm c;
  comm_init(&c, comm);

  struct crystal cr;
  crystal_init(&cr, &c);

  sarray_transfer(elem_data, &elements, proc, 0, &cr);
  *nelt_ = nelt = elements.n;
  ed = elements.ptr;

  vl = *vl_ = (long long *)realloc(*vl_, nv * nelt * sizeof(long long));
  for (e = 0; e < nelt; ++e)
    for (n = 0; n < nv; ++n)
      vl[e * nv + n] = ed[e].vtx[n];

  if (coord != NULL) {
    coord = *coord_ =
        (double *)realloc(*coord_, ndim * nv * nelt * sizeof(double));
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

void parrsb_get_part_stat(int *nc, int *ns, int *nss, int *nel, long long *vtx,
                          int nelt, int nv, MPI_Comm ce) {
  struct comm comm;
  comm_init(&comm, ce);

  int np = comm.np;
  int id = comm.id;

  if (np == 1)
    return;

  int Npts = nelt * nv;
  int i;
  slong *data = (slong *)malloc((Npts + 1) * sizeof(slong));
  for (i = 0; i < Npts; i++)
    data[i] = vtx[i];
  struct gs_data *gsh = gs_setup(data, Npts, &comm, 0, gs_pairwise, 0);

  int Nmsg;
  pw_data_nmsg(gsh, &Nmsg);

  int *Ncomm = (int *)malloc((Nmsg + 1) * sizeof(int));
  pw_data_size(gsh, Ncomm);

  gs_free(gsh);
  free(data);

  int nelMin, nelMax, nelSum;
  int ncMin, ncMax, ncSum;
  int nsMin, nsMax, nsSum;
  int nssMin, nssMax, nssSum;
  int b;

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

  nelMax = nelt;
  nelMin = nelt;
  nelSum = nelt;
  comm_allreduce(&comm, gs_int, gs_max, &nelMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nelMin, 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &nelSum, 1, &b);

  free(Ncomm);
  comm_free(&comm);

  if (nc != NULL) {
    nc[0] = ncMin;
    nc[1] = ncMax;
    nc[2] = ncSum;
  }

  if (ns != NULL) {
    ns[0] = nsMin;
    ns[1] = nsMax;
    ns[2] = nsSum;
  }

  if (nss != NULL) {
    nss[0] = nssMin;
    nss[1] = nssMax;
    nss[2] = nssSum;
  }

  if (nel != NULL) {
    nel[0] = nelMin;
    nel[1] = nelMax;
    nel[2] = nelSum;
  }
}

void parrsb_print_part_stat(long long *vtx, int nelt, int nv, MPI_Comm ce) {
  int id, np;
  MPI_Comm_rank(ce, &id);
  MPI_Comm_size(ce, &np);

  int nc[3], ns[3], nss[3], nel[3];
  parrsb_get_part_stat(&nc[0], &ns[0], &nss[0], &nel[0], vtx, nelt, nv, ce);

  if (id == 0) {
    printf("Min neighbors: %d | Max neighbors: %d | Avg neighbors: %lf\n",
           nc[0], nc[1], (double)nc[2] / np);
    printf("Min nvolume: %d | Max nvolume: %d | Avg nvolume: %lf\n", ns[0],
           ns[1], (double)ns[2] / np);
    printf("Min volume: %d | Max volume: %d | Avg volume: %lf\n", nss[0],
           nss[1], (double)nss[2] / np);
    printf("Min elements: %d | Max elements: %d\n", nel[0], nel[1]);
    fflush(stdout);
  }
}

struct parrsb_input *parrsb_parse_input(int argc, char *argv[], MPI_Comm comm) {
  struct parrsb_input *in = tcalloc(struct parrsb_input, 1);
  in->mesh = NULL, in->tol = 0.2, in->test = 0, in->dump = 1;
  in->verbose = 0, in->type = 0;
  MPI_Comm_size(comm, &in->nactive);

  static struct option long_options[] = {{"mesh", required_argument, 0, 'm'},
                                         {"tol", optional_argument, 0, 't'},
                                         {"test", no_argument, 0, 'c'},
                                         {"no-dump", no_argument, 0, 'd'},
                                         {"nactive", optional_argument, 0, 'n'},
                                         {"verbose", optional_argument, 0, 'v'},
                                         {"type", optional_argument, 0, 'y'},
                                         {0, 0, 0, 0}};

  for (;;) {
    int opt_idx = 0, c = getopt_long(argc, argv, "", long_options, &opt_idx);
    if (c == -1)
      break;

    switch (c) {
    case 'm':
      in->mesh = optarg;
      break;
    case 't':
      in->tol = atof(optarg);
      break;
    case 'c':
      in->test = 1;
      break;
    case 'd':
      in->dump = 0;
      break;
    case 'n':
      in->nactive = atoi(optarg);
      break;
    case 'v':
      in->verbose = atoi(optarg);
      break;
    case 'y':
      in->type = atoi(optarg);
      break;
    case '?':
      break;
    default:
      exit(1);
    }
  }

  return in;
}

void parrsb_check_error_(int err, char *file, int line, MPI_Comm comm) {
  int sum;
  MPI_Allreduce(&err, &sum, 1, MPI_INT, MPI_SUM, comm);

  if (sum != 0) {
    int id;
    MPI_Comm_rank(comm, &id);
    if (id == 0)
      printf("check_error failure in %s:%d\n", file, line);

    MPI_Finalize();
    exit(1);
  }
}
