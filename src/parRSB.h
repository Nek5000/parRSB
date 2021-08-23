#ifndef _PARRSB_H_
#define _PARRSB_H_

#include "gslib.h"

#if !defined(MPI)
#error "gslib needs to be compiled with MPI"
#endif

#if !defined(GLOBAL_LONG_LONG)
#error "gslib needs to be compiled with GLOBAL_LONG_LONG"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Partitioning
 */
typedef struct {
  /* General options */
  int partitioner;   // 0 - RSB, 1 - RCB, 2 - RIB (Default: 0)
  int debug_level;   // 0, 1, 2, .. etc (Default: 0)
  int profile_level; // 0, 1, 2, .. etc (Default: 0)

  /* RSB specific */
  int rsb_algo;     // 0 - Lanczos, 1 - RQI (Default: 0)
  int rsb_pre;      // 0 - None, 1 - RCB , 2 - RIB (Default: 1)
  int rsb_grammian; // 0 or 1 (Default: 1)

  /* Other */
  int repair; // 0 - No, 1 - Yes (Default: 1)
} parrsb_options;

extern parrsb_options parrsb_default_options;

int parRSB_partMesh(int *part, int *seq, long long *vtx, double *coord, int nel,
                    int nv, parrsb_options options, MPI_Comm comm);

#define fparRSB_partMesh FORTRAN_UNPREFIXED(fparrsb_partmesh, FPARRSB_PARTMESH)
void fparRSB_partMesh(int *part, int *seq, long long *vtx, double *coord,
                      int *nel, int *nve, int *options, int *comm, int *err);

/*
 * Connectivity
 */
int parRSB_findConnectivity(long long *vtx, double *coord, int nel, int nDim,
                            long long *periodicInfo, int nPeriodicFaces,
                            double tol, MPI_Comm comm, int verbose);

#define fparRSB_findConnectivity                                               \
  FORTRAN_UNPREFIXED(fparrsb_findconnectivity, FPARRSB_FINDCONNECTIVITY)
void fparRSB_findConnectivity(long long *vtx, double *coord, int *nel,
                              int *nDim, long long *periodicInfo,
                              int *nPeriodicFaces, double *tol, MPI_Fint *fcomm,
                              int *verbose, int *err);

/*
 * Auxiliary functions
 */
typedef struct {
  char *mesh; // Mesh name, required.
  double tol; // gencon tolerance, default: 0.2
  int test; // run tests, default: 0
  int dump; // dump the connectivity or map file, default: 1
  int nactive; // # of active MPI ranks, default: MPI_Comm_size
} parrsb_input;

int parrsb_read_mesh(unsigned int *nel, int *nv, long long **vl, double **coord,
                     unsigned int *nbcs, long long **bcs, char *name,
                     MPI_Comm comm, int read);

int parrsb_dump_con(long long *vl, unsigned int nelt, int nv, char *name,
                    MPI_Comm comm);

int parrsb_dump_map(int nelt, int nv, int *pmap, long long *vtx, char *name,
                    MPI_Comm comm);

int parrsb_distribute_elements(unsigned int *nelt, long long **vl,
                               double **coord, int *part, int nv,
                               MPI_Comm comm);

void parrsb_print_part_stat(long long *vtx, int nelt, int nv, MPI_Comm comm);

void parrsb_get_part_stat(int *nc, int *ns, int *nss, int *nel, long long *vtx,
                          int nelt, int nv, MPI_Comm comm);

parrsb_input *parrsb_parse_input(int argc, char *argv[]);

void parrsb_check_error_(int err, char *file, int line, MPI_Comm comm);
#define parrsb_check_error(err, comm) parrsb_check_error_(err, __FILE__, __LINE__, comm)

#ifdef __cplusplus
}
#endif

#endif
