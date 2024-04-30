#ifndef _PARRSB_METRICS_H_
#define _PARRSB_METRICS_H_

#include "gslib.h"

//------------------------------------------------------------------------------
// Metrics
//
typedef enum {
  RSB_BALANCE = 0,
  RSB_COMPONENTS,
  RSB_COMPONENTS_NCOMP,
  RSB_FIEDLER,
  RSB_FIEDLER_SETUP,
  RSB_FIEDLER_CALC,
  RSB_FIEDLER_CALC_NITER,
  RSB_INVERSE_SETUP,
  RSB_INVERSE,
  RSB_LANCZOS_SETUP,
  RSB_LANCZOS,
  RSB_LANCZOS_TQLI,
  RSB_NEIGHBORS,
  RSB_PRE,
  RSB_PROJECT_AX,
  RSB_PROJECT_MG,
  RSB_SORT,
  TOL_FNL,
  TOL_TGT,
  TOL_INIT
} metric;

void metric_init(void);
void metric_acc(metric m, double val);
void metric_set(metric m, double val);
void metric_tic(struct comm *c, metric m);
void metric_toc(struct comm *c, metric m);
double metric_get_value(int level, metric m);
void metric_push_level(void);
uint metric_get_levels(void);
void metric_rsb_print(struct comm *c, int profile_level);
void metric_crs_print(struct comm *c, int profile_level);
void metric_finalize(void);

#endif
