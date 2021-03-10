#ifndef _GENMAP_H_
#define _GENMAP_H_

#include "genmap-gslib.h"
#include "genmap-types.h"
#include <mpi.h>

#include <parRSB.h>

typedef struct comm *genmap_comm;
typedef struct genmap_handle_private *genmap_handle;
typedef struct genmap_vector_private *genmap_vector;
typedef struct rsb_element *genmap_element;

int genmap_init(genmap_handle *h, comm_ext ce, parRSB_options *options);

genmap_element genmap_get_elements(genmap_handle h);

genmap_comm genmap_local_comm(genmap_handle h);
genmap_comm genmap_global_comm(genmap_handle h);

void genmap_set_nvertices(genmap_handle h, int nVertices);
void genmap_set_elements(genmap_handle h, struct array *localElements);

GenmapULong genmap_get_partition_nel(genmap_handle h);
void genmap_set_partition_nel(genmap_handle h, GenmapULong globalElements);

GenmapLong genmap_get_local_start_index(genmap_handle h);
void genmap_set_local_start_index(genmap_handle h, GenmapLong localStart);

GenmapInt genmap_get_nel(genmap_handle h);

int genmap_finalize(genmap_handle h);

/* genmap_comm */
void genmap_comm_scan(genmap_handle h, genmap_comm c);
int genmap_comm_size(genmap_comm c);
int genmap_comm_rank(genmap_comm c);
void genmap_comm_split(struct comm *old, int bin, int key, struct comm *new_);

/* Genamp Vector */
int GenmapCreateVector(genmap_vector *x, GenmapInt size);
int GenmapSetVector(genmap_vector x, GenmapScalar *array);
int GenmapGetVector(genmap_vector x, GenmapScalar *array);

int GenmapCreateRandomVector(genmap_vector *x, GenmapInt size, GenmapInt seed);
int GenmapCreateOnesVector(genmap_vector *x, GenmapInt size);
int GenmapCreateZerosVector(genmap_vector *x, GenmapInt size);

int GenmapScaleVector(genmap_vector y, genmap_vector x, GenmapScalar alpha);
int GenmapAxpbyVector(genmap_vector z, genmap_vector x, GenmapScalar alpha,
                      genmap_vector y, GenmapScalar beta);

int GenmapCopyVector(genmap_vector x, genmap_vector y);
GenmapScalar GenmapDotVector(genmap_vector x, genmap_vector y);

int GenmapOrthogonalizebyOneVector(struct comm *c, genmap_vector q1,
                                   GenmapULong n);

int GenmapPrintVector(genmap_vector x);
int GenmapDestroyVector(genmap_vector x);

/* Laplacian */
int GenmapInitLaplacian(genmap_handle h, struct comm *c);
int GenmapLaplacian(genmap_handle h, GenmapScalar *u, GenmapScalar *v);

int GenmapInitLaplacianWeighted(genmap_handle h, struct comm *c);
int GenmapLaplacianWeighted(genmap_handle h, GenmapScalar *u, GenmapScalar *v);

/* Eigen */
int GenmapTQLI(genmap_handle h, genmap_vector diag, genmap_vector upper,
               genmap_vector **eVec, genmap_vector *eVal);
int genmap_inverse_power(double *y, int N, double *A, int verbose);
int genmap_power(double *y, int N, double *A, int verbose);

/* Matrix inverse */
void matrix_inverse(int N, double *A);

/* Lanczos */
int GenmapLanczosLegendary(genmap_handle h, struct comm *c, genmap_vector f,
                           GenmapInt niter, genmap_vector **rr,
                           genmap_vector diag, genmap_vector upper);
int GenmapLanczos(genmap_handle h, struct comm *c, genmap_vector init,
                  GenmapInt iter, genmap_vector **q, genmap_vector alpha,
                  genmap_vector beta);

/* Fiedler */
int GenmapFiedlerLanczos(genmap_handle h, struct comm *c, int maxIter,
                         int global);
int GenmapFiedlerRQI(genmap_handle h, struct comm *c, int maxIter, int global);

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
