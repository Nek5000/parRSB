#ifndef _GENMAP_PRECON_H_
#define _GENMAP_PRECON_H_

#include <genmap-impl.h>

struct mg_data *mg_setup(int factor, struct comm *c, struct csr_laplacian *M);
void mg_vcycle(GenmapScalar *u, GenmapScalar *rhs, struct mg_data *d);
void mg_free(struct mg_data *d);

#endif
