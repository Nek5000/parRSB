#include "metrics.h"
#include "gslib.h"
#include <limits.h>
#include <time.h>

#define MAXMETS 50
#define MAXLVLS 30
#define MAXSIZE (MAXMETS * MAXLVLS)

static double metrics[MAXMETS];
static double *stack;
static uint stack_size;

void metric_init(void) {
  for (uint i = 0; i < MAXMETS; i++) metrics[i] = 0.0;
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
  if (level < 0) return metrics[m];
  if ((uint)level < stack_size) return stack[level * MAXMETS + m];
  return 0.0;
}

void metric_push_level(void) {
  assert(stack_size < MAXLVLS && "stack_size >= MAXLVLS");

  for (unsigned i = 0; i < MAXMETS; i++) {
    stack[stack_size * MAXMETS + i] = metrics[i];
    metrics[i] = 0.0;
  }
  stack_size++;
}

uint metric_get_levels(void) { return stack_size; }

static void metric_print_aux(double *wrk, struct comm *c) {
  double *min = wrk, *max = min + MAXSIZE, *sum = max + MAXSIZE;
  double *buf = sum + MAXSIZE;

  uint max_size = stack_size * MAXMETS;
  for (uint i = 0; i < max_size; i++) { min[i] = max[i] = sum[i] = stack[i]; }

  comm_allreduce(c, gs_double, gs_min, min, MAXSIZE, buf); // min
  comm_allreduce(c, gs_double, gs_max, max, MAXSIZE, buf); // max
  comm_allreduce(c, gs_double, gs_add, sum, MAXSIZE, buf); // sum
  for (uint i = 0; i < max_size; i++) sum[i] /= c->np;
}

#define SUMMARY(i, m)                                                          \
  sum[i * MAXMETS + m], min[i * MAXMETS + m], max[i * MAXMETS + m]

void metric_rsb_print(struct comm *c, int profile_level) {
  double *wrk = tcalloc(double, 4 * MAXSIZE);
  metric_print_aux(wrk, c);
  double *min = wrk, *max = min + MAXSIZE, *sum = max + MAXSIZE;

  uint i;
  for (i = 0; i < stack_size; i++) {
    if (c->id == 0 && profile_level > 0) {
      printf("level=%02d\n", i);
      printf("  RSB_PRE                    : %e/%e/%e\n", SUMMARY(i, RSB_PRE));
      printf("  RSB_FIEDLER                : %e/%e/%e\n",
             SUMMARY(i, RSB_FIEDLER));
      printf("    RSB_FIEDLER_SETUP        : %e/%e/%e\n",
             SUMMARY(i, RSB_FIEDLER_SETUP));
      printf("    RSB_FIEDLER_CALC         : %e/%e/%e\n",
             SUMMARY(i, RSB_FIEDLER_CALC));
      printf("      RSB_LANCZOS_SETUP      : %e/%e/%e\n",
             SUMMARY(i, RSB_LANCZOS_SETUP));
      printf("      RSB_LANCZOS            : %e/%e/%e\n",
             SUMMARY(i, RSB_LANCZOS));
      printf("      RSB_LANCZOS_TQLI       : %e/%e/%e\n",
             SUMMARY(i, RSB_LANCZOS_TQLI));
      printf("    RSB_FIEDLER_CALC_NITER   : %e/%e/%e\n",
             SUMMARY(i, RSB_FIEDLER_CALC_NITER));
      printf("  RSB_SORT                   : %e/%e/%e\n", SUMMARY(i, RSB_SORT));
      printf("  RSB_COMPONENTS             : %e/%e/%e\n",
             SUMMARY(i, RSB_COMPONENTS));
      printf("    RSB_COMPONENTS_NCOMP     : %e/%e/%e\n",
             SUMMARY(i, RSB_COMPONENTS_NCOMP));
      printf("  RSB_NEIGHBORS              : %e/%e/%e\n",
             SUMMARY(i, RSB_NEIGHBORS));
      printf("  RSB_BALANCE                : %e/%e/%e\n",
             SUMMARY(i, RSB_BALANCE));
    }
    fflush(stdout);
  }

  if (wrk) free(wrk);
}

#undef SUMMARY

void metric_finalize(void) {
  if (stack != NULL) free(stack), stack = NULL;
}

#undef MAXMETS
#undef MAXLVLS
#undef MAXSIZE
