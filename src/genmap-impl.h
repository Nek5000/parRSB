#ifndef _GENMAP_IMPL_H_
#define _GENMAP_IMPL_H_

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include <genmap-multigrid-precon.h>
#include <genmap.h>

#define GENMAP_ALIGN 32

#define GENMAP_SP_TOL 1e-05
#define GENMAP_DP_TOL 1e-12
#define GENMAP_TOL GENMAP_DP_TOL

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

int rcb(struct array *elements, size_t unit_size, int ndim, struct comm *c,
        buffer *bfr);
int rib(struct array *elements, size_t unit_size, int ndim, struct comm *c,
        buffer *bfr);

/* Laplacian */
struct gs_laplacian {
  struct gs_data *gsh;
  GenmapScalar *weights;
  GenmapScalar *u;
  uint lelt;
  int nv;
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

/* Fiedler */
int GenmapFiedlerLanczos(struct rsb_element *elements, uint lelt, int nv,
                         int max_iter, int global, struct comm *gsc,
                         buffer *buf);

int GenmapFiedlerRQI(struct rsb_element *elements, uint lelt, int nv,
                     int max_iter, int global, struct comm *gsc, buffer *buf);

/* RSB */
int genmap_rsb(genmap_handle h);

/* Iterative methods */
int flex_cg(genmap_vector x, struct gs_laplacian *gl, mgData d,
            genmap_vector ri, int maxIter, struct comm *gsc, buffer *buf);

int project(genmap_vector x, struct gs_laplacian *gl, mgData d,
            genmap_vector ri, int max_iter, struct comm *gsc, buffer *buf);

struct genmap_handle_private {
  struct comm *global;

  GenmapLong nel;
  GenmapLong Nnodes;
  GenmapLong start;
  int nv;

  struct array *elements;
  struct crystal cr;

  /* Weighted Laplacian */
  struct gs_data *gsw;
  GenmapScalar *weights;
  buffer buf;

  parrsb_options *options;
  size_t unit_size;
};

struct genmap_vector_private {
  GenmapInt size;
  GenmapScalar *data;
};

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

/* Repair and balance */
sint get_components(sint *component, struct rsb_element *elements,
                    struct comm *c, buffer *buf, uint nelt, uint nv);

int repair_partitions(genmap_handle h, struct comm *tc, struct comm *lc,
                      int bin, struct comm *gc);
int balance_partitions(genmap_handle h, struct comm *lc, int bin,
                       struct comm *gc);

/* Matrix inverse */
void matrix_inverse(int N, double *A);

/* Dump data */
int GenmapFiedlerDump(const char *fname, genmap_handle h, struct comm *c);
int GenmapVectorDump(const char *fname, GenmapScalar *y, genmap_handle h,
                     struct comm *c);
int GenmapCentroidDump(const char *fname, genmap_handle h, sint g_id,
                       struct comm *c);
int GenmapElementIdDump(const char *fname, genmap_handle h, struct comm *c);

#endif
