#include "metrics.h"
#include "gslib.h"
#include <limits.h>
#include <time.h>

#define MAXMETS 5100
#define MAXLVLS 20
#define MAXSIZE (MAXMETS * MAXLVLS)

static double metrics[MAXMETS];
static double *stack;
static uint stack_size;

void metric_init() {
  for (uint i = 0; i < MAXMETS; i++)
    metrics[i] = 0.0;
  stack = tcalloc(double, MAXSIZE);
  stack_size = 0;
}

void metric_acc(metric m, double val) { metrics[m] += val; }

void metric_set(metric m, double val) { metrics[m] = val; }

void metric_tic(struct comm *c, metric m) {
  comm_barrier(c);
  metrics[m] -= comm_time();
}

void metric_toc(struct comm *c, metric m) {
  metrics[m] += comm_time();
  comm_barrier(c);
}

double metric_get_value(int level, metric m) {
  if (level == -1)
    return metrics[m];
  else if (level < stack_size)
    return stack[level * MAXMETS + m];
  return 0.0;
}

void metric_push_level() {
  assert(stack_size < MAXLVLS && "stack_size >= MAXLVLS");

  uint i;
  for (i = 0; i < MAXMETS; i++) {
    stack[stack_size * MAXMETS + i] = metrics[i];
    metrics[i] = 0.0;
  }
  stack_size++;
}

uint metric_get_levels() { return stack_size; }

static void metric_print_aux(double **mini, struct comm *c, int profile_level) {
  *mini = tcalloc(double, 4 * MAXSIZE);
  double *min = *mini, *max = min + MAXSIZE, *sum = max + MAXSIZE,
         *buf = sum + MAXSIZE;

  uint max_size = stack_size * MAXMETS, i;
  for (i = 0; i < max_size; i++)
    min[i] = max[i] = sum[i] = stack[i];

  comm_allreduce(c, gs_double, gs_min, min, MAXSIZE, buf); // min
  comm_allreduce(c, gs_double, gs_max, max, MAXSIZE, buf); // max
  comm_allreduce(c, gs_double, gs_add, sum, MAXSIZE, buf); // sum
  for (i = 0; i < max_size; i++)
    sum[i] /= c->np;
}

#define SUMMARY(i, m)                                                          \
  sum[i * MAXMETS + m], min[i * MAXMETS + m], max[i * MAXMETS + m]

void metric_rsb_print(struct comm *c, int profile_level) {
  double *min = NULL;
  metric_print_aux(&min, c, profile_level);
  double *max = min + MAXSIZE, *sum = max + MAXSIZE;

  uint i;
  for (i = 0; i < stack_size; i++) {
    if (c->id == 0 && profile_level > 0) {
      printf("level=%02d\n", i);
      printf("  RSB_PRE                    : %e/%e/%e\n", SUMMARY(i, RSB_PRE));
      printf("  RSB_FIEDLER                : %e/%e/%e\n",
             SUMMARY(i, RSB_FIEDLER));
      printf("  RSB_FIEDLER_SETUP          : %e/%e/%e\n",
             SUMMARY(i, RSB_FIEDLER_SETUP));
      printf("  RSB_FIEDLER_INNER          : %e/%e/%e\n",
             SUMMARY(i, RSB_FIEDLER_INNER));
      printf("    RSB_LANCZOS_SETUP        : %e/%e/%e\n",
             SUMMARY(i, RSB_LANCZOS_SETUP));
      printf("    RSB_LANCZOS              : %e/%e/%e\n",
             SUMMARY(i, RSB_LANCZOS));
      printf("    RSB_LANCZOS_TQLI         : %e/%e/%e\n",
             SUMMARY(i, RSB_LANCZOS_TQLI));
      printf("    RSB_INVERSE_SETUP        : %e/%e/%e\n",
             SUMMARY(i, RSB_INVERSE_SETUP));
      printf("    RSB_INVERSE              : %e/%e/%e\n",
             SUMMARY(i, RSB_INVERSE));
      printf("  RSB_FIEDLER_INNER_NITER    : %e/%e/%e\n",
             SUMMARY(i, RSB_FIEDLER_INNER_NITER));
      printf("  RSB_SORT                   : %e/%e/%e\n", SUMMARY(i, RSB_SORT));
      printf("  RSB_REPAIR_BALANCE         : %e/%e/%e\n",
             SUMMARY(i, RSB_REPAIR_BALANCE));
    }
  }

  if (min)
    free(min);
}

void metric_crs_print(struct comm *c, int profile_level) {
  double *min = NULL;
  metric_print_aux(&min, c, profile_level);
  double *max = min + MAXSIZE, *sum = max + MAXSIZE;

  uint i;
  for (i = 0; i < stack_size; i++) {
    if (c->id == 0 && profile_level > 0) {
      printf("level=%02d\n", i);
      printf("  SCHUR_SOLVE_CHOL1  : %e/%e/%e\n",
             SUMMARY(i, SCHUR_SOLVE_CHOL1));
      printf("  SCHUR_SOLVE_SETRHS1: %e/%e/%e\n",
             SUMMARY(i, SCHUR_SOLVE_SETRHS1));
      printf("  SCHUR_SOLVE_PROJECT: %e/%e/%e\n",
             SUMMARY(i, SCHUR_SOLVE_PROJECT));
      printf("\tSCHUR_PROJECT_NITER: %e/%e/%e\n",
             SUMMARY(i, SCHUR_PROJECT_NITER));
      printf("\tSCHUR_PROJECT_OPERATOR: %e/%e/%e\n",
             SUMMARY(i, SCHUR_PROJECT_OPERATOR));
      printf("\tSCHUR_PROJECT_PRECOND : %e/%e/%e\n",
             SUMMARY(i, SCHUR_PROJECT_PRECOND));
      printf("  SCHUR_SOLVE_SETRHS2: %e/%e/%e\n",
             SUMMARY(i, SCHUR_SOLVE_SETRHS2));
      printf("  SCHUR_SOLVE_CHOL2  : %e/%e/%e\n",
             SUMMARY(i, SCHUR_SOLVE_CHOL2));
    }
  }

  if (min)
    free(min);
}

#undef SUMMARY

void metric_finalize() {
  if (stack != NULL)
    free(stack), stack = NULL;
}

#undef MAXMETS
#undef MAXLVLS
#undef MAXSIZE
