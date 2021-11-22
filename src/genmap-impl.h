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

#define MAXDIM 3 /* Maximum dimension of the mesh */
#define MAXNV 8  /* Maximum number of vertices per element */

/* `struct rcb_element` is used for RCB and RIB partitioning.
   `struct rsb_element` should be a superset of `struct rcb_element` */
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

/* CSR Matrix */
struct csr_mat {
  uint rn;
  ulong row_start;

  uint *roff;
  ulong *col;
  GenmapScalar *v;
  GenmapScalar *diag;
  GenmapScalar *buf;

  struct gs_data *gsh;
};

void csr_mat_setup(struct csr_mat *M, struct array *entries, struct comm *c,
                   buffer *buf);
void csr_mat_gather(GenmapScalar *buf, struct csr_mat *M, GenmapScalar *x,
                    buffer *bfr);
void csr_mat_apply(GenmapScalar *y, struct csr_mat *M, GenmapScalar *x,
                   buffer *buf);
void csr_mat_print(struct csr_mat *M, struct comm *c);
int csr_mat_free(struct csr_mat *M);

/* RSB */
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

#define GS 1
#define CSR 2
#define WEIGHTED 4
#define UNWEIGHTED 8

struct vector {
  GenmapInt size;
  GenmapScalar *data;
};
typedef struct vector *genmap_vector;

struct laplacian {
  int type;
  uint nel;
  int nv;

  /* GS */
  GenmapScalar *diag;
  struct gs_data *gsh;

  /* CSR */
  struct csr_mat *M;

  /* GPU */
  /* unique column ids of local laplacian matrix, sorted */
  ulong *col_ids;
  /* adj as csr, values are not stored as everything is -1 */
  uint *adj_off;
  uint *adj_ind;
  /* diagonal as an array */
  uint *diag_ind;
  uint *diag_val;

  /* WORK ARRAY
   * size = nel * nv for gs, # of unique column ids for csr
   */
  GenmapScalar *u;
};

int laplacian_init(struct laplacian *l, struct rsb_element *elems, uint nel,
                   int nv, int type, struct comm *c, buffer *buf);
int laplacian(GenmapScalar *v, struct laplacian *l, GenmapScalar *u,
              buffer *buf);
void laplacian_free(struct laplacian *l);

void matrix_inverse(int N, double *A);

int power_serial(double *y, int N, double *A, int verbose);

int fiedler(struct rsb_element *elements, uint nel, int nv, int max_iter,
            int global, struct comm *gsc, buffer *buf, int gid);

int repair_partitions(struct array *elements, int nv, struct comm *tc,
                      struct comm *lc, int bin, struct comm *gc, buffer *bfr);
int balance_partitions(struct array *elements, int nv, struct comm *lc,
                       struct comm *gc, int bin, buffer *bfr);

int rsb(struct array *elements, parrsb_options *options, int nv,
        struct comm *gc, buffer *bfr);

/* occa */

int occa_init(char *backend, int device_id, int platform_id);
int occa_lanczos_init(struct comm *c, struct laplacian *l, int niter);
int occa_lanczos(GenmapScalar *diag, GenmapScalar *upper, GenmapScalar *rr,
                 GenmapScalar *f, struct comm *c, buffer *bfr);
int occa_lanczos_free();

typedef struct {
  GenmapULong sequenceId;
  int nNeighbors;
  GenmapULong elementId;
  GenmapULong vertexId;
  uint workProc;
} vertex;

/* Memory */
int GenmapMallocArray(size_t n, size_t unit, void *p);
int GenmapCallocArray(size_t n, size_t unit, void *p);
int GenmapReallocArray(size_t n, size_t unit, void *p);
int GenmapFree(void *p);

#define GenmapMalloc(n, p) GenmapMallocArray((n), sizeof(**(p)), p)
#define GenmapCalloc(n, p) GenmapCallocArray((n), sizeof(**(p)), p)
#define GenmapRealloc(n, p) GenmapReallocArray((n), sizeof(**(p)), p)

/* Genmap Metrics */
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

/* Dump data */
int GenmapFiedlerDump(const char *fname, struct rsb_element *elm, uint nelt,
                      int nv, struct comm *c);
int GenmapVectorDump(const char *fname, GenmapScalar *y,
                     struct rsb_element *elm, uint nelt, int nv,
                     struct comm *c);
int GenmapElementDump(const char *fname, struct rsb_element *elm, uint nelt,
                      int nv, struct comm *c, int dump);
/* Misc */
int log2ll(long long n);
int logbll(long long n, int a);

void genmap_barrier(struct comm *c);
void comm_split(struct comm *old, int bin, int key, struct comm *new_);

double GenmapGetMaxRss();
void GenmapPrintStack();

#endif
