/*
 * Generate connectivity (.co2) from Nek5000 mesh file (.re2).
 */
#include <getopt.h>
#include <mpi.h>
#include <stdio.h>

#include <parRSB.h>

struct input {
  char *mesh;
  double tol;
  int test;
  int dump;
  int nactive;
};

static int parse_input(struct input *in, int argc, char *argv[]) {
  in->mesh = NULL;
  in->tol = 0.2;
  in->test = 0;
  in->dump = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &in->nactive);

  int c;
  for (;;) {
    static struct option long_options[] = {
        {"mesh", required_argument, 0, 'm'},
        {"tol", optional_argument, 0, 't'},
        {"test", no_argument, 0, 'c'},
        {"no-dump", no_argument, 0, 'd'},
        {"nactive", optional_argument, 0, 'n'},
        {0, 0, 0, 0}};

    int option_index = 0;
    c = getopt_long(argc, argv, "", long_options, &option_index);

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
    case '?':
      break;
    default:
      abort();
    }
  }
}

static void check_error_(int err, char *file, int line, MPI_Comm comm) {
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

#define check_error(err) check_error_(err, __FILE__, __LINE__, MPI_COMM_WORLD);

struct pair {
  ulong gid;
  uint indx;
};

static int test_parcon(unsigned int neltp, long long *vlp, char *name,
                       MPI_Comm comm) {
  unsigned int nelt;
  int nv;
  long long *vls = NULL;
  int err = parrsb_read_mesh(&nelt, &nv, &vls, NULL, NULL, NULL, name,
                             MPI_COMM_WORLD, 2);
  assert(neltp == nelt);

  uint size = nelt * nv;
  slong *minp = tcalloc(slong, size);
  slong *maxp = tcalloc(slong, size);

  uint i;
  for (i = 0; i < size; i++)
    minp[i] = maxp[i] = vlp[i];

  buffer buf;
  buffer_init(&buf, 1024);

  struct comm c;
  comm_init(&c, comm);

  struct gs_data *gsh = gs_setup(vls, size, &c, 0, gs_pairwise, 0);
  gs(minp, gs_long, gs_min, 0, gsh, &buf);
  gs(maxp, gs_long, gs_max, 0, gsh, &buf);

  for (i = 0; i < size; i++)
    if (minp[i] != maxp[i]) {
      err = 1;
      break;
    }

  int np;
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  if (np == 1) {
    for (i = 0; i < size; i++)
      if (vls[i] != vlp[i]) {
        err = 1;
        break;
      }
  }

  /* TODO: Add a check for number of distinct global ids serial vs parallel and
   * that should be enough to make sure serial vs parallel mapping the same upto
   * a permutation. Or set the gs using vlp and check min/max of vls */

  gs_free(gsh);
  comm_free(&c);
  buffer_free(&buf);

  if (minp != NULL)
    free(minp);
  if (maxp != NULL)
    free(maxp);
  if (vls != NULL)
    free(vls);

  return err;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  struct input in;
  parse_input(&in, argc, argv);

  /* Read the geometry from the .re2 file */
  unsigned int nelt, nbcs;
  double *coord = NULL;
  long long *bcs = NULL;
  int nv;
  int err = parrsb_read_mesh(&nelt, &nv, NULL, &coord, &nbcs, &bcs, in.mesh,
                             MPI_COMM_WORLD, 1);
  check_error(err);

  /* Find connectivity */
  long long *vl = (long long *)calloc(nelt * nv, sizeof(long long));
  int ndim = nv == 8 ? 3 : 2;
  err |= parRSB_findConnectivity(vl, coord, nelt, ndim, bcs, nbcs, in.tol,
                                 MPI_COMM_WORLD, 0);
  check_error(err);

  /* Test if the relevant env. variable is set */
  if (in.test == 1)
    err |= test_parcon(nelt, vl, in.mesh, MPI_COMM_WORLD);
  check_error(err);

  /* Write connectivity file */
  if (in.dump == 1)
    err |= parrsb_dump_con(vl, nelt, nv, in.mesh, MPI_COMM_WORLD);
  check_error(err);

  /* Free resources */
  if (vl != NULL)
    free(vl);
  if (coord != NULL)
    free(coord);
  if (bcs != NULL)
    free(bcs);

  MPI_Finalize();

  return 0;
}

#undef check_error
