#ifndef _GENMAP_H_
#define _GENMAP_H_

#include <genmap-types.h>
#include <parRSB.h>

typedef struct genmap_handle_private *genmap_handle;
typedef struct genmap_vector_private *genmap_vector;

/* genmap_handle */
int genmap_init(genmap_handle *h, comm_ext ce, parrsb_options *options);
void *genmap_get_elements(genmap_handle h);
int genmap_get_nvertices(genmap_handle h);
GenmapULong genmap_get_partition_nel(genmap_handle h);
GenmapInt genmap_get_nel(genmap_handle h);
GenmapLong genmap_get_local_start_index(genmap_handle h);
int genmap_finalize(genmap_handle h);

void genmap_comm_scan(genmap_handle h, struct comm *c);
void genmap_barrier(struct comm *c);
void genmap_comm_split(struct comm *old, int bin, int key, struct comm *new_);

/* genmap_vector */
int genmap_vector_create(genmap_vector *x, GenmapInt size);
int genmap_vector_create_zeros(genmap_vector *x, GenmapInt size);
int genmap_vector_copy(genmap_vector x, genmap_vector y);

int genmap_vector_scale(genmap_vector y, genmap_vector x, GenmapScalar alpha);
int genmap_vector_axpby(genmap_vector z, genmap_vector x, GenmapScalar alpha,
                        genmap_vector y, GenmapScalar beta);

GenmapScalar genmap_vector_dot(genmap_vector x, genmap_vector y);

int genmap_vector_ortho_one(struct comm *c, genmap_vector q1, GenmapULong n);

int genmap_destroy_vector(genmap_vector x);

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

/* Lanczos */
int GenmapLanczos(genmap_handle h, struct comm *c, genmap_vector f, int niter,
                  genmap_vector **rr, genmap_vector diag, genmap_vector upper);

/* Fiedler */
int GenmapFiedlerLanczos(genmap_handle h, struct comm *c, int maxIter,
                         int global);
int GenmapFiedlerRQI(genmap_handle h, struct comm *c, int maxIter, int global);

/* RSB/RCB */
int genmap_rsb(genmap_handle h);

/* Misc */
double GenmapGetMaxRss();
void GenmapPrintStack();
int log2ll(long long n);

#endif
