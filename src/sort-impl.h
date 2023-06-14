#ifndef _PARRSB_SORT_IMPL_H_
#define _PARRSB_SORT_IMPL_H_

#include "sort.h"

double get_scalar(struct array *a, uint i, uint offset, uint usize,
                  gs_dom type);

void get_extrema(void *extrema_, struct sort *data, uint field,
                 const struct comm *c);

uint *set_proc_from_idx(uint size, sint np, slong start, slong nelem);

void sort_local(struct sort *s);

void sarray_transfer_chunk(struct array *arr, const size_t usize,
                           const uint *proc, const struct comm *c);

struct hypercube {
  struct sort *data;
  int nprobes;
  double *probes;
  ulong *probe_cnt;
};
void parallel_hypercube_sort(struct hypercube *data, const struct comm *c);

void parallel_bin_sort(struct sort *s, const struct comm *c);

#endif // _PARRSB_SORT_IMPL_H_
