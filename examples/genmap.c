/*
 * Generate partitions (.ma2) from Nek5000's mesh (.re2) file.
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

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  struct input in;
  parse_input(&in, argc, argv);

  int id;
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  int color = 0;
  if (id < in.nactive)
    color = 1;
  MPI_Comm comm;
  MPI_Comm_split(MPI_COMM_WORLD, color, id, &comm);

  /* Read the geometry from the .re2 file */
  unsigned int nelt, nbcs;
  double *coord = NULL;
  long long *bcs = NULL;
  int nv;
  int err = 0;
  if (color == 1)
    err = parrsb_read_mesh(&nelt, &nv, NULL, &coord, &nbcs, &bcs, in.mesh, comm,
                           1);
  check_error(err);

  /* Find connectivity */
  long long *vl = (long long *)calloc(nelt * nv, sizeof(long long));
  int ndim = nv == 8 ? 3 : 2;
  if (color == 1)
    err |= parRSB_findConnectivity(vl, coord, nelt, ndim, bcs, nbcs, in.tol,
                                   comm, 0);
  check_error(err);

  if (color == 1)
    parrsb_print_part_stat(vl, nelt, nv, comm);

  /* Partition the mesh */
  parrsb_options options = parrsb_default_options;
  int *part = (int *)calloc(nelt, sizeof(int));
  if (color == 1)
    err |= parRSB_partMesh(part, NULL, vl, coord, nelt, nv, options, comm);
  check_error(err);

  /* Redistribute data */
  if (color == 1)
    err |= parrsb_distribute_elements(&nelt, &vl, &coord, part, nv, comm);
  check_error(err);

  if (color == 1)
    parrsb_print_part_stat(vl, nelt, nv, comm);

  /* Write map file */
  if (color == 1 && in.dump == 1)
    err |= parrsb_dump_map(nelt, nv, part, vl, in.mesh, comm);
  check_error(err);

  /* Free resources */
  if (part != NULL)
    free(part);
  if (vl != NULL)
    free(vl);
  if (coord != NULL)
    free(coord);
  if (bcs != NULL)
    free(bcs);

  MPI_Comm_free(&comm);
  MPI_Finalize();

  return 0;
}

#undef check_error
