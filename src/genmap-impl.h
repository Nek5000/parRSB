#ifndef _GENMAP_IMPL_H_
#define _GENMAP_IMPL_H_

#include "parRSB.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef scalar
#undef scalar
#endif
#define scalar double

#ifdef SCALAR_MAX
#undef SCALAR_MAX
#endif
#define SCALAR_MAX DBL_MAX

#ifdef gs_scalar
#undef gs_scalar
#endif
#define gs_scalar gs_double

#define GENMAP_ALIGN 32
#define GENMAP_TOL 1e-12

#define MAXDIM 3 // Maximum dimension of the mesh
#define MAXNV 8  // Maximum number of vertices per element

//------------------------------------------------------------------------------
// RCB / RIB
// `struct rcb_element` is used for RCB and RIB partitioning.
// `struct rsb_element` should be a superset of `struct rcb_element`
struct rcb_element {
  uint proc, origin, seq;
  ulong globalId;
  scalar coord[MAXDIM], fiedler;
};

int rib(struct array *elements, size_t unit_size, int ndim, struct comm *c,
        buffer *bfr);

//------------------------------------------------------------------------------
// RSB
//
struct rsb_element {
  uint proc, origin, seq;
  ulong globalId;
  scalar coord[MAXDIM], fiedler;
  slong vertices[MAXNV];
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
int laplacian(scalar *v, struct laplacian *l, scalar *u, buffer *buf);
void laplacian_free(struct laplacian *l);

//------------------------------------------------------------------------------
// OCCA
int occa_init(const char *backend, int device_id, int platform_id,
              struct comm *c);
int occa_lanczos_init(struct comm *c, struct laplacian *l, int niter);
int occa_lanczos_aux(scalar *diag, scalar *upper, scalar *rr, uint lelt,
                     ulong nelg, int niter, scalar *f, struct laplacian *gl,
                     struct comm *gsc, buffer *bfr);
int occa_lanczos_free();
int occa_free();

//------------------------------------------------------------------------------
// Misc
//
int log2ll(long long n);
void genmap_barrier(struct comm *c);

#endif
