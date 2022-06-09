#ifndef _GENMAP_METRICS_H_
#define _GENMAP_METRICS_H_

#include "gslib.h"

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
  SCHUR_PROJECT_OPERATOR,
  SCHUR_PROJECT_PRECOND,
  SCHUR_SOLVE_CHOL1,
  SCHUR_SOLVE_CHOL2,
  SCHUR_SOLVE_PROJECT,
  SCHUR_SOLVE_SETRHS1,
  SCHUR_SOLVE_SETRHS2,
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
void metric_rsb_print(struct comm *c, int profile_level);
void metric_crs_print(struct comm *c, int profile_level);
void metric_finalize();

#endif
