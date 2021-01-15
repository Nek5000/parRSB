#ifndef _GENMAP_IMPL_H_
#define _GENMAP_IMPL_H_

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#ifdef GENMAP_DEBUG
#include <stdio.h>
#endif

#include <genmap-multigrid-precon.h>
#include <genmap.h>

#define GENMAP_FIEDLER 0
#define GENMAP_GLOBALID 1
#define GENMAP_PROC 2
#define GENMAP_ORIGIN 3

#define GENMAP_RCB_ELEMENT 0
#define GENMAP_RSB_ELEMENT 1

#define GENMAP_SUM 0
#define GENMAP_MAX 1
#define GENMAP_MIN 2
#define GENMAP_MUL 3

#define GENMAP_ALIGN 32

#define GENMAP_SP_TOL 1e-05
#define GENMAP_DP_TOL 1e-12
#define GENMAP_TOL GENMAP_DP_TOL

#define GENMAP_READER_LEN 256
#define GENMAP_MAX_READERS 32

#define MAXDIM 3 /* Maximum dimension of the mesh */
#define MAXNV 8  /* Maximum number of vertices per element */

struct genmap_comm_private {
  struct comm gsc;

  /* Un-weighted Laplacian */
  struct gs_data *gsh;
  csr_mat M;

  /* Weighted Laplacian */
  struct gs_data *gsw;
  buffer buf;
  GenmapScalar *b;
};

/* parRCB/parRSB internals */
struct rcb_element {
  unsigned char type;
  GenmapInt proc;
  GenmapInt origin;
  GenmapInt seq;
  GenmapLong globalId;
  GenmapScalar coord[MAXDIM];
  GenmapScalar fiedler;
};

int rcb(struct comm *ci, struct array *elements, int ndim, buffer *bfr);
int rib(struct comm *ci, struct array *elements, int ndim, buffer *bfr);

/* rsb_element should be a superset of rcb_element */
struct rsb_element {
  unsigned char type;
  GenmapInt proc;
  GenmapInt origin;
  GenmapInt seq;
  GenmapLong globalId;
  GenmapScalar coord[MAXDIM];
  GenmapScalar fiedler;
  GenmapLong vertices[8];
  GenmapInt part;
  GenmapULong globalId0;
};

int GenmapCreateElements(GenmapElements *e);
int GenmapDestroyElements(GenmapElements e);
GenmapElements GenmapGetElements_default(genmap_handle h);

struct genmap_handle_private {
  genmap_comm global;
  genmap_comm local;

  GenmapLong nel;
  GenmapLong Nnodes;
  GenmapLong start;
  int nv;

  GenmapScalar *weights;
  struct array *elements;
  struct crystal cr;

  parRSB_options *options;
};

int GenmapCreateHandle(genmap_handle h);
int GenmapDestroyHandle(genmap_handle h);

struct GenmapVector_private {
  GenmapInt size;
  GenmapScalar *data;
};

#define GenmapMalloc(n, p) GenmapMallocArray((n), sizeof(**(p)), p)
#define GenmapCalloc(n, p) GenmapCallocArray((n), sizeof(**(p)), p)
#define GenmapRealloc(n, p) GenmapReallocArray((n), sizeof(**(p)), p)

/* Genmap Metrics */
typedef enum {
  RCB,
  WEIGHTEDLAPLACIANSETUP,
  FIEDLER,
  NFIEDLER,
  FIEDLERSORT,
  BISECT,
  LANCZOS,
  NLANCZOS,
  WEIGHTEDLAPLACIAN,
  TQLI,
  LAPLACIANSETUP,
  PRECONDSETUP,
  RQI,
  NRQI,
  PROJECT,
  NPROJECT,
  GRAMMIAN,
  LAPLACIAN,
  VCYCLE,
  END
} metric;

void metric_init();
void metric_finalize();
void metric_acc(metric m, double count);
void metric_tic(struct comm *c, metric m);
void metric_toc(struct comm *c, metric m);
double metric_get_value(int level, metric m);
void metric_push_level();
uint metric_get_levels();
void metric_print(struct comm *c);

/* genCon */
typedef struct {
  GenmapULong sequenceId;
  int nNeighbors;
  GenmapULong elementId;
  GenmapULong vertexId;
  uint workProc;
} vertex;

/* Components */
sint is_disconnected(struct comm *c, struct gs_data *gsh, buffer *buf,
                     uint nelt, uint nv);

/* Matrix inverse */
void matrix_inverse(int N, double *A);

/* Dump data */
int GenmapFiedlerDump(const char *fname, genmap_handle h, struct comm *c);
int GenmapVectorDump(const char *fname, GenmapScalar *y, uint size,
                     struct comm *c);
int GenmapCentroidDump(const char *fname, genmap_handle h, struct comm *c);

#endif
