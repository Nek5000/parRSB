#include <stdio.h>
#include <stdlib.h>

#include <genmap-impl.h>

void genmap_barrier(struct comm *c) {
#if defined(GENMAP_SYNC_BY_REDUCTION)
  sint dummy = c->id, wrk;
  comm_allreduce(c, gs_int, gs_max, &dummy, 1, &wrk);
#else
  comm_barrier(c);
#endif
}

int log2ll(long long n) {
  int k = 0;
  while (n > 1)
    n /= 2, k++;

  return k;
}
