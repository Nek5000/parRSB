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

#define GENMAP_SP_TOL 1e-05
#define GENMAP_DP_TOL 1e-12
#define GENMAP_TOL GENMAP_DP_TOL

#define MAXDIM 3 /* Maximum dimension of the mesh */
#define MAXNV 8  /* Maximum number of vertices per element */

/* genmap_vector */
struct genmap_vector_private {
  GenmapInt size;
  GenmapScalar *data;
};

typedef struct genmap_vector_private *genmap_vector;

int genmap_vector_create(genmap_vector *x, GenmapInt size);
int genmap_vector_create_zeros(genmap_vector *x, GenmapInt size);
int genmap_vector_copy(genmap_vector x, genmap_vector y);

int genmap_vector_scale(genmap_vector y, genmap_vector x, GenmapScalar alpha);
int genmap_vector_axpby(genmap_vector z, genmap_vector x, GenmapScalar alpha,
                        genmap_vector y, GenmapScalar beta);

GenmapScalar genmap_vector_dot(genmap_vector x, genmap_vector y);
int genmap_vector_ortho_one(struct comm *c, genmap_vector q1, GenmapULong n);
int genmap_destroy_vector(genmap_vector x);

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

struct rsb_element {
  GenmapInt proc;
  GenmapInt origin;
  GenmapInt seq;
  GenmapULong globalId;
  GenmapScalar coord[MAXDIM];
  GenmapScalar fiedler;
  GenmapLong vertices[8];
  GenmapInt part;
};

/* Laplacian */
struct gs_laplacian {
  struct gs_data *gsh;
  GenmapScalar *weights;
  GenmapScalar *u;
  uint lelt;
  int nv;
};

typedef struct {
  ulong r, c;
  uint proc;
} csr_entry;

typedef struct {
  ulong r, c, rn, cn;
  uint p;
  GenmapScalar v;
} entry;

struct csr_mat {
  uint rn;
  ulong row_start;

  uint *roff;
  ulong *col;
  GenmapScalar *v, *diag, *buf;

  struct gs_data *gsh;
};

int GenmapInitLaplacian(struct csr_mat *M, struct rsb_element *elems, uint n,
                        int nv, struct comm *c, buffer *buf);
int GenmapLaplacian(GenmapScalar *v, struct csr_mat *M, GenmapScalar *u,
                    buffer *buf);

int GenmapInitLaplacianWeighted(struct gs_laplacian *gl,
                                struct rsb_element *elems, uint lelt, int nv,
                                struct comm *c, buffer *buf);
int GenmapLaplacianWeighted(GenmapScalar *v, struct gs_laplacian *gl,
                            GenmapScalar *u, buffer *buf);

/* Eigen */
int genmap_inverse_power(double *y, int N, double *A, int verbose);
int genmap_power(double *y, int N, double *A, int verbose);

/* Fiedler */
int GenmapFiedlerLanczos(struct rsb_element *elements, uint lelt, int nv,
                         int max_iter, int global, struct comm *gsc,
                         buffer *buf);

int GenmapFiedlerRQI(struct rsb_element *elements, uint lelt, int nv,
                     int max_iter, int global, struct comm *gsc, buffer *buf);

/* Repair and balance */
int repair_partitions(struct array *elements, int nv, struct comm *tc,
                      struct comm *lc, int bin, struct comm *gc, buffer *bfr);
int balance_partitions(struct array *elements, int nv, struct comm *lc,
                       struct comm *gc, int bin, buffer *bfr);

/* RSB */
int rsb(struct array *elements, parrsb_options *options, int nv,
        struct comm *gc, buffer *bfr);

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
  RCB,
  FIEDLER,
  FIEDLER_NITER,
  LANCZOS,
  LANCZOS_NITER,
  LANCZOS_TOL_FINAL,
  LANCZOS_TOL_TARGET,
  LAPLACIAN,
  LAPLACIAN_INIT,
  RQI,
  RQI_NITER,
  PROJECT,
  PROJECT_NITER,
  COMPONENTS,
  END
} metric;

void metric_init();
void metric_acc(metric m, double count);
void metric_tic(struct comm *c, metric m);
void metric_toc(struct comm *c, metric m);
double metric_get_value(int level, metric m);
void metric_push_level();
uint metric_get_levels();
void metric_print(struct comm *c, int profile_level);
void metric_finalize();

/* Gencon */
typedef struct {
  GenmapULong sequenceId;
  int nNeighbors;
  GenmapULong elementId;
  GenmapULong vertexId;
  uint workProc;
} vertex;

/* Matrix inverse */
void matrix_inverse(int N, double *A);

/* Dump data */
int GenmapFiedlerDump(const char *fname, struct rsb_element *elm, uint nelt,
                      int nv, struct comm *c);
int GenmapVectorDump(const char *fname, GenmapScalar *y,
                     struct rsb_element *elm, uint nelt, int nv,
                     struct comm *c);
/* Misc */
void genmap_barrier(struct comm *c);
void genmap_comm_split(struct comm *old, int bin, int key, struct comm *new_);
double GenmapGetMaxRss();
void GenmapPrintStack();
int log2ll(long long n);

#endif
