#ifndef _GENMAP_IMPL_H_
#define _GENMAP_IMPL_H_

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include <genmap-types.h>
#include <parRSB.h>

#define GENMAP_ALIGN 32
#define GENMAP_TOL 1e-12

#define MAXDIM 3 // Maximum dimension of the mesh
#define MAXNV 8  // Maximum number of vertices per element

//------------------------------------------------------------------------------
// RCB / RIB
// `struct rcb_element` is used for RCB and RIB partitioning.
// `struct rsb_element` should be a superset of `struct rcb_element`
struct rcb_element {
  GenmapInt proc;
  GenmapInt origin;
  GenmapInt seq;
  GenmapULong globalId;
  GenmapScalar coord[MAXDIM];
  GenmapScalar fiedler;
};

int rcb(struct array *elements, size_t unit_size, int ndim, struct comm *c,
        buffer *bfr);
int rib(struct array *elements, size_t unit_size, int ndim, struct comm *c,
        buffer *bfr);

//------------------------------------------------------------------------------
// RSB
//
struct rsb_element {
  GenmapInt proc;
  GenmapInt origin;
  GenmapInt seq;
  GenmapULong globalId;
  GenmapScalar coord[MAXDIM];
  GenmapScalar fiedler;
  GenmapLong vertices[MAXNV];
  GenmapInt part;
};

//------------------------------------------------------------------------------
// Laplacian
//
#define GS 1
#define CSR 2
#define GPU 4
#define WEIGHTED 1024
#define UNWEIGHTED 2048

struct laplacian {
  int type, nv;
  uint nel;
  void *data;
};

struct csr_laplacian {
  uint rn;
  ulong rstart;
  uint *roff;

  ulong *col;

  GenmapScalar *v;

  GenmapScalar *diag;

  GenmapScalar *buf;

  struct gs_data *gsh;
};

void csr_mat_gather(GenmapScalar *buf, struct csr_laplacian *M, GenmapScalar *x,
                    buffer *bfr);
void csr_mat_apply(GenmapScalar *y, struct csr_laplacian *M, GenmapScalar *x,
                   buffer *buf);
int csr_mat_free(struct csr_laplacian *M);

struct gs_laplacian {
  GenmapScalar *diag, *u;
  struct gs_data *gsh;
};

struct gpu_laplacian {
  uint rn;

  // unique column ids of local laplacian matrix, sorted
  ulong *col_ids;

  // adj as csr, for unweighted case, adj_val is null as the values are all -1
  uint *adj_off;
  uint *adj_ind;
  GenmapScalar *adj_val;

  /* diagonal as an array */
  uint *diag_ind;
  GenmapScalar *diag_val;

  // gs for host side communication
  struct gs_data *gsh;
};

int laplacian_init(struct laplacian *l, struct rsb_element *elems, uint nel,
                   int nv, int type, struct comm *c, buffer *buf);
int laplacian(GenmapScalar *v, struct laplacian *l, GenmapScalar *u,
              buffer *buf);
void laplacian_free(struct laplacian *l);

//------------------------------------------------------------------------------
// Repair and balance
int repair_partitions(struct array *elements, int nv, struct comm *tc,
                      struct comm *lc, int bin, struct comm *gc, buffer *bfr);
int balance_partitions(struct array *elements, int nv, struct comm *lc,
                       struct comm *gc, int bin, buffer *bfr);

//------------------------------------------------------------------------------
// RSB
int fiedler(struct rsb_element *elements, uint nel, int nv, int max_iter,
            int global, struct comm *gsc, buffer *buf);

int rsb(struct array *elements, parrsb_options *options, int nv,
        struct comm *gc, buffer *bfr);

//------------------------------------------------------------------------------
// OCCA
typedef struct vector *genmap_vector;
int occa_init(char *backend, int device_id, int platform_id);
int occa_lanczos_init(struct comm *c, struct laplacian *l, int niter);
int occa_lanczos_aux(genmap_vector diag, genmap_vector upper, genmap_vector *rr,
                     uint lelt, ulong nelg, int niter, genmap_vector f,
                     struct laplacian *gl, struct comm *gsc, buffer *bfr);
int occa_lanczos_free();
int occa_free();

//------------------------------------------------------------------------------
// Memory
//
int GenmapMallocArray(size_t n, size_t unit, void *p);
int GenmapCallocArray(size_t n, size_t unit, void *p);
int GenmapReallocArray(size_t n, size_t unit, void *p);
int GenmapFree(void *p);

#define GenmapMalloc(n, p) GenmapMallocArray((n), sizeof(**(p)), p)
#define GenmapCalloc(n, p) GenmapCallocArray((n), sizeof(**(p)), p)
#define GenmapRealloc(n, p) GenmapReallocArray((n), sizeof(**(p)), p)

//------------------------------------------------------------------------------
// Metrics
//
typedef enum {
  PRE,
  FIEDLER,
  FIEDLER_NITER,
  FIEDLER_SORT,
  REPAIR_BALANCE,
  LANCZOS,
  TOL_FINAL,
  TOL_TARGET,
  LAPLACIAN,
  LAPLACIAN_INIT,
  RQI,
  PROJECT,
  PROJECT_NITER,
  COMPONENTS,
  END
} metric;

void metric_init();
void metric_acc(metric m, double val);
void metric_set(metric m, double val);
void metric_tic(struct comm *c, metric m);
void metric_toc(struct comm *c, metric m);
double metric_get_value(int level, metric m);
void metric_push_level();
uint metric_get_levels();
void metric_print(struct comm *c, int profile_level);
void metric_finalize();

//------------------------------------------------------------------------------
// Misc
//
struct vector {
  GenmapInt size;
  GenmapScalar *data;
};

typedef struct {
  GenmapULong sequenceId;
  int nNeighbors;
  GenmapULong elementId;
  GenmapULong vertexId;
  uint workProc;
} vertex;

void matrix_inverse(int N, double *A);
int power_serial(double *y, int N, double *A, int verbose);

int GenmapFiedlerDump(const char *fname, struct rsb_element *elm, uint nelt,
                      int nv, struct comm *c);
int GenmapVectorDump(const char *fname, GenmapScalar *y,
                     struct rsb_element *elm, uint nelt, int nv,
                     struct comm *c);
int GenmapElementDump(const char *fname, struct rsb_element *elm, uint nelt,
                      int nv, struct comm *c, int dump);
int log2ll(long long n);
int logbll(long long n, int a);

void genmap_barrier(struct comm *c);
void comm_split(struct comm *old, int bin, int key, struct comm *new_);

double GenmapGetMaxRss();
void GenmapPrintStack();

#endif
