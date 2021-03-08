#ifndef _GENMAP_H_
#define _GENMAP_H_

#include "genmap-gslib.h"
#include "genmap-types.h"
#include <mpi.h>

#include <parRSB.h>

typedef MPI_Datatype GenmapDataType;
typedef struct genmap_comm_private *genmap_comm;
typedef struct genmap_handle_private *genmap_handle;
typedef struct GenmapVector_private *GenmapVector;
typedef struct rsb_element *GenmapElements;

int genmap_init(genmap_handle *h, comm_ext ce, parRSB_options *options);
int genmap_finalize(genmap_handle h);

int GenmapMallocArray(size_t n, size_t unit, void *p);
int GenmapCallocArray(size_t n, size_t unit, void *p);
int GenmapReallocArray(size_t n, size_t unit, void *p);
int GenmapFree(void *p);

GenmapElements GenmapGetElements(genmap_handle h);

genmap_comm GenmapGetLocalComm(genmap_handle h);
void GenmapSetLocalComm(genmap_handle h, genmap_comm c);

genmap_comm GenmapGetGlobalComm(genmap_handle h);
void GenmapSetGlobalComm(genmap_handle h, genmap_comm c);

GenmapInt GenmapGetNLocalElements(genmap_handle h);
void genmap_set_elements(genmap_handle h, struct array *localElements);

GenmapLong genmap_get_global_nel(genmap_handle h);
void GenmapSetNGlobalElements(genmap_handle h, GenmapLong globalElements);

GenmapLong GenmapGetLocalStartIndex(genmap_handle h);
void GenmapSetLocalStartIndex(genmap_handle h, GenmapLong localStart);

int GenmapGetNVertices(genmap_handle h);
void genmap_set_vertices(genmap_handle h, int nVertices);

void genmap_scan(genmap_handle h, genmap_comm c);
void genmap_comm_scan(genmap_handle h, struct comm *c);

int GenmapCreateComm(genmap_comm *c, comm_ext ce);
int GenmapDestroyComm(genmap_comm c);

int genmap_comm_size(genmap_comm c);
int genmap_comm_rank(genmap_comm c);

int GenmapGop(genmap_comm c, void *v, GenmapInt size, GenmapDataType type,
              GenmapInt op);
int GenmapReduce(genmap_comm c, void *out, void *in, GenmapInt size,
                 GenmapDataType type, GenmapInt op);
int GenmapBcast(genmap_comm c, void *in, GenmapInt count, GenmapDataType type);

void comm_split(struct comm *old, int bin, int key, struct comm *new_);
void GenmapSplitComm(genmap_handle h, genmap_comm *c, int bin);

int GenmapCrystalInit(genmap_handle h, genmap_comm c);
int GenmapCrystalTransfer(genmap_handle h, int field);
int GenmapCrystalFinalize(genmap_handle h);

int GenmapRead(genmap_handle h, void *data);

/* Genamp Vector */
int GenmapCreateVector(GenmapVector *x, GenmapInt size);
int GenmapSetVector(GenmapVector x, GenmapScalar *array);
int GenmapGetVector(GenmapVector x, GenmapScalar *array);

int GenmapCreateRandomVector(GenmapVector *x, GenmapInt size, GenmapInt seed);
int GenmapCreateOnesVector(GenmapVector *x, GenmapInt size);
int GenmapCreateZerosVector(GenmapVector *x, GenmapInt size);

int GenmapScaleVector(GenmapVector y, GenmapVector x, GenmapScalar alpha);
int GenmapAxpbyVector(GenmapVector z, GenmapVector x, GenmapScalar alpha,
                      GenmapVector y, GenmapScalar beta);

int GenmapVectorsEqual(GenmapVector x, GenmapVector y, GenmapScalar tol);
int GenmapCopyVector(GenmapVector x, GenmapVector y);
GenmapScalar GenmapDotVector(GenmapVector x, GenmapVector y);
GenmapScalar GenmapAbsMaxVector(GenmapVector x);
GenmapScalar GenmapMaxVector(GenmapVector x);
GenmapScalar GenmapAbsMinVector(GenmapVector x);
GenmapScalar GenmapMinVector(GenmapVector x);
GenmapScalar GenmapNormVector(GenmapVector x, GenmapInt p);

int GenmapOrthogonalizebyOneVector(genmap_comm c, GenmapVector q1,
                                   GenmapLong n);

int GenmapPrintVector(GenmapVector x);
int GenmapDestroyVector(GenmapVector x);

/* Laplacian */
int GenmapInitLaplacian(genmap_handle h, struct comm *c);
int GenmapLaplacian(genmap_handle h, GenmapScalar *u, GenmapScalar *v);

int GenmapInitLaplacianWeighted(genmap_handle h, struct comm *c);
int GenmapLaplacianWeighted(genmap_handle h, GenmapScalar *u, GenmapScalar *v);

/* Eigen */
int GenmapTQLI(genmap_handle h, GenmapVector diag, GenmapVector upper,
               GenmapVector **eVec, GenmapVector *eVal);
int genmap_inverse_power(double *y, int N, double *A, int verbose);
int genmap_power(double *y, int N, double *A, int verbose);

/* Matrix inverse */
void matrix_inverse(int N, double *A);

/* Lanczos */
int GenmapLanczosLegendary(genmap_handle h, genmap_comm c, GenmapVector f,
                           GenmapInt niter, GenmapVector **rr,
                           GenmapVector diag, GenmapVector upper);
int GenmapLanczos(genmap_handle h, genmap_comm c, GenmapVector init,
                  GenmapInt iter, GenmapVector **q, GenmapVector alpha,
                  GenmapVector beta);

/* Fiedler */
int GenmapFiedlerLanczos(genmap_handle h, genmap_comm c, int maxIter,
                         int global);
int GenmapFiedlerRQI(genmap_handle h, genmap_comm c, int maxIter, int global);

/* RSB/RCB */
void genmap_load_balance(struct array *eList, uint nel, int nv, double *coord,
                         long long *vtx, struct crystal *cr, buffer *bfr);
int genmap_rsb(genmap_handle h);
int genmap_rcb(genmap_handle h);
int genmap_rib(genmap_handle h);

/* Misc */
double GenmapGetMaxRss();
void GenmapPrintStack();

#endif
