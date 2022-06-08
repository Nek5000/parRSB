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
  GenmapInt proc, origin, seq;
  GenmapULong globalId;
  GenmapScalar coord[MAXDIM], fiedler;
};

int rcb(struct array *elements, size_t unit_size, int ndim, struct comm *c,
        buffer *bfr);
int rib(struct array *elements, size_t unit_size, int ndim, struct comm *c,
        buffer *bfr);

//------------------------------------------------------------------------------
// RSB
//
struct rsb_element {
  GenmapInt proc, origin, seq;
  GenmapULong globalId;
  GenmapScalar coord[MAXDIM], fiedler;
  GenmapLong vertices[MAXNV];
  GenmapInt part;
};

int rsb(struct array *elements, parrsb_options *options, int nv,
        struct comm *gc, buffer *bfr);

//------------------------------------------------------------------------------
// Laplacian
//
#define GS 1
#define CSR 2
#define CSC 4

struct laplacian;
struct laplacian *laplacian_init(struct rsb_element *elems, uint nel, int nv,
                                 int type, struct comm *c, buffer *buf);
int laplacian(GenmapScalar *v, struct laplacian *l, GenmapScalar *u,
              buffer *buf);
void laplacian_free(struct laplacian *l);

//------------------------------------------------------------------------------
// OCCA
int occa_init(const char *backend, int device_id, int platform_id,
              struct comm *c);
int occa_lanczos_init(struct comm *c, struct laplacian *l, int niter);
int occa_lanczos_aux(GenmapScalar *diag, GenmapScalar *upper, GenmapScalar *rr,
                     uint lelt, ulong nelg, int niter, GenmapScalar *f,
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
  RSB_COMPONENTS,
  RSB_FIEDLER,
  RSB_FIEDLER_NITER,
  RSB_FIEDLER_SORT,
  RSB_LANCZOS,
  RSB_LAPLACIAN,
  RSB_LAPLACIAN_INIT,
  RSB_PRE,
  RSB_PROJECT,
  RSB_PROJECT_NITER,
  RSB_REPAIR_BALANCE,
  TOL_FNL,
  TOL_TGT
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
int log2ll(long long n);
void genmap_barrier(struct comm *c);
void comm_split(const struct comm *old, int bin, int key, struct comm *new_);

#endif
