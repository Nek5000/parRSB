#ifndef _GENMAP_IMPL_H_
#define _GENMAP_IMPL_H_

#include "genmap-types.h"
#include "parRSB.h"
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

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

int rib(struct array *elements, size_t unit_size, int ndim, struct comm *c,
        buffer *bfr);

//------------------------------------------------------------------------------
// RSB
//
struct rsb_element {
  GenmapInt proc, origin, seq;
  GenmapULong globalId;
  GenmapScalar coord[MAXDIM], fiedler;
  slong vertices[MAXNV];
  GenmapInt part;
};

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
// Misc
//
int log2ll(long long n);
void genmap_barrier(struct comm *c);
void comm_split(const struct comm *old, int bin, int key, struct comm *new_);

#endif
