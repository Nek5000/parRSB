#include <float.h>
#include "sort.h"

#define min(a, b) ((a) < (b) ? (a) : (b))

double get_scalar(struct array *a, uint i, uint offset, uint usize,
                  gs_dom type) {
  char *v = (char *)a->ptr + i * usize + offset;

  double data;
  switch (type) {
  case gs_int:
    data = *((uint *)v);
    break;
  case gs_long:
    data = *((ulong *)v);
    break;
  case gs_double:
    data = *((double *)v);
    break;
  default:
    break;
  }

  return data;
}

void get_extrema(void *extrema_, struct sort *data, uint field,
                 struct comm *c) {
  struct array *a = data->a;
  uint usize = data->unit_size;
  uint offset = data->offset[field];
  gs_dom t = data->t[field];

  double *extrema = extrema_;
  sint size = a->n;
  if (size == 0) {
    extrema[0] = -DBL_MAX;
    extrema[1] = -DBL_MAX;
  } else {
    extrema[0] = -get_scalar(a, 0, offset, usize, t);
    extrema[1] = get_scalar(a, size - 1, offset, usize, t);
  }

  double buf[2];
  comm_allreduce(c, gs_double, gs_max, extrema, 2, buf);
  extrema[0] *= -1;
}

int set_dest(uint *proc, uint size, sint np, slong start, slong nelem) {
  uint i;
  if (nelem < np) {
    for (i = 0; i < size; i++)
      proc[i] = start + i;
    return 0;
  }

  uint psize = nelem / np;
  uint nrem = nelem - np * psize;
  slong lower = nrem * (psize + 1);

  uint id1, id2;
  if (start < lower)
    id1 = start / (psize + 1);
  else
    id1 = nrem + (start - lower) / psize;

  if (start + size < lower)
    id2 = (start + size) / (psize + 1);
  else
    id2 = nrem + (start + size - lower) / psize;

  i = 0;
  while (id1 <= id2 && i < size) {
    ulong s = id1 * psize + min(id1, nrem);
    ulong e = (id1 + 1) * psize + min(id1 + 1, nrem);
    e = min(e, nelem);
    assert(id1 < np && "Ooops ! id1 is larger than np");
    assert(id1 >= 0 && "Ooops ! id1 is smaller than 0");
    while (s <= start + i && start + i < e && i < size)
      proc[i++] = id1;
    id1++;
  }

  return 0;
}

//-----------------------------------------------------------------------------
// Parallel Bin-Sort
//
int set_bin(uint **proc_, struct sort *s, uint field, struct comm *c) {
  struct array *a = s->a;
  gs_dom t = s->t[field];
  uint offset = s->offset[field];

  uint size = a->n;
  uint *proc = *proc_ = tcalloc(uint, size);

  double extrema[2];
  get_extrema((void *)extrema, s, field, c);
  double range = extrema[1] - extrema[0];

  if (size == 0)
    return 0;

  sint np = c->np;
  uint id = 0;
  uint index = 0;
  do {
    double end = extrema[0] + (range / np) * (id + 1);
    while (index < size) {
      double val = get_scalar(a, index, offset, s->unit_size, t);
      if (val <= end)
        proc[index++] = id;
      else
        break;
    }
    id++;
  } while (id < np && index < size);
  for (; index < size; index++)
    proc[index] = np - 1;
}

int parallel_bin_sort(struct sort *s, struct comm *c) {
  // Local sort
  sort_local(s);

  // Set destination bin
  uint *proc;
  set_bin(&proc, s, 0, c);

  // Transfer to destination processor
  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer_ext_(s->a, s->unit_size, proc, sizeof(uint), &cr);
  crystal_free(&cr);

  free(proc);

  // Locally sort again
  sort_local(s);

  return 0;
}

//-----------------------------------------------------------------------------
// Parallel Hypercube-Sort
//
struct hypercube {
  struct sort *data;
  int nprobes;
  double *probes;
  ulong *probe_cnt;
};

int init_probes(struct hypercube *data, struct comm *c) {
  struct sort *input = data->data;

  // Allocate space for probes and counts
  int nprobes = data->nprobes = 3;
  if (!data->probes)
    data->probes = tcalloc(double, nprobes);
  if (!data->probe_cnt)
    data->probe_cnt = tcalloc(ulong, nprobes);

  double extrema[2];
  get_extrema((void *)extrema, data->data, 0, c);
  double range = extrema[1] - extrema[0];
  double delta = range / (nprobes - 1);

  data->probes[0] = extrema[0];
  data->probes[1] = extrema[0] + delta;
  data->probes[2] = extrema[1];

  return 0;
}

int update_probe_counts(struct hypercube *data, struct comm *c) {
  struct sort *input = data->data;
  uint offset = input->offset[0];
  gs_dom t = input->t[0];

  uint nprobes = data->nprobes;
  uint i;
  for (i = 0; i < nprobes; i++)
    data->probe_cnt[i] = 0;

  struct array *a = input->a;
  uint e;
  for (e = 0; e < a->n; e++) {
    double val_e = get_scalar(a, e, offset, input->unit_size, t);
    for (i = 0; i < nprobes; i++)
      if (val_e < data->probes[i])
        data->probe_cnt[i]++;
  }

  ulong buf[3];
  comm_allreduce(c, gs_long, gs_add, data->probe_cnt, nprobes, buf);

  return 0;
}

int update_probes(slong nelem, double *probes, ulong *probe_cnt,
                  uint threshold) {
  slong expected = nelem / 2;
  if (llabs(expected - probe_cnt[1]) < threshold)
    return 0;

  if (probe_cnt[1] < expected)
    probes[0] = probes[1];
  else
    probes[2] = probes[1];

  probes[1] = probes[0] + (probes[2] - probes[0]) / 2.0;

  return 0;
}

int transfer_elem(struct hypercube *data, struct comm *c) {
  struct sort *input = data->data;
  uint usize = input->unit_size;
  uint offset = input->offset[0];
  gs_dom t = input->t[0];

  struct array *a = input->a;
  uint size = a->n;

  uint e, lown = 0, uppern = 0;
  for (e = 0; e < size; e++) {
    double val = get_scalar(a, e, offset, usize, t);
    if (val < data->probes[1])
      lown++;
    else
      uppern++;
  }

  slong out[2][2], in[2], buf[2][2];
  in[0] = lown, in[1] = uppern;
  comm_scan(out, c, gs_long, gs_add, in, 2, buf);

  ulong lstart = out[0][0], ustart = out[0][1];
  ulong lelem = out[1][0], uelem = out[1][1];
  uint np = c->np, lnp = np / 2, unp = np - lnp;

  uint *proc = tcalloc(uint, size);
  set_dest(proc, lnp, lstart, lown, lelem);
  set_dest(proc + lown, unp, ustart, uppern, uelem);

  for (e = lown; e < size; e++)
    proc[e] += lnp;

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer_ext_(a, usize, proc, sizeof(uint), &cr);
  crystal_free(&cr);

  free(proc);

  return 0;
}

int parallel_hypercube_sort(struct hypercube *data, struct comm *c) {
  struct sort *input = data->data;
  struct array *a = input->a;
  gs_dom t = input->t[0];
  uint offset = input->offset[0];

  sint size = c->np, rank = c->id;

  slong out[2][1], buf[2][1], in = a->n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
  slong start = out[0][0];
  slong nelem = out[1][0];

  uint threshold = nelem / (10 * size);
  if (threshold < 2)
    threshold = 2;

  sort_local(data->data);

  if (size == 1)
    return 0;

  init_probes(data, c);
  update_probe_counts(data, c);

  int max_iter = log2((data->probes[2] - data->probes[0]) / 1e-12);
  int iter = 0;
  while (llabs(nelem / 2 - data->probe_cnt[1]) > threshold &&
         iter++ < max_iter) {
    update_probes(nelem, data->probes, data->probe_cnt, threshold);
    update_probe_counts(data, c);
  }

  transfer_elem(data, c);

  // split the communicator
  struct comm nc;
  sint lower = (rank < size / 2) ? 1 : 0;
#if defined(MPI)
  MPI_Comm nc_;
  MPI_Comm_split(c->c, lower, rank, &nc_);
  comm_init(&nc, nc_);
  MPI_Comm_free(&nc_);
#else
  comm_init(&nc, 1);
#endif

  // TODO: Keep load balancing after each split
  parallel_hypercube_sort(data, &nc);
  comm_free(&nc);

  return 0;
}

int load_balance(struct array *a, size_t size, struct comm *c,
                 struct crystal *cr) {
  slong out[2][1], buf[2][1];
  slong in = a->n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
  slong start = out[0][0];
  slong nelem = out[1][0];

  uint *proc = tcalloc(uint, a->n);
  set_dest(proc, a->n, c->np, start, nelem);
  sarray_transfer_ext_(a, size, proc, sizeof(uint), cr);
  free(proc);

  return 0;
}

int parallel_sort_private(struct sort *data, struct comm *c) {
  struct comm dup;
  comm_dup(&dup, c);

  int balance = data->balance;
  sort_algo algo = data->algo;

  struct array *a = data->a;
  size_t usize = data->unit_size;

  struct hypercube hdata;

  switch (algo) {
  case bin_sort:
    parallel_bin_sort(data, c);
    break;
  case hypercube_sort:
    hdata.data = data;
    hdata.probes = NULL;
    hdata.probe_cnt = NULL;
    parallel_hypercube_sort(&hdata, &dup);
    free(hdata.probes);
    free(hdata.probe_cnt);
    break;
  default:
    break;
  }

  if (balance) {
    struct crystal cr;
    crystal_init(&cr, c);
    load_balance(a, usize, c, &cr);
    crystal_free(&cr);
    sort_local(data);
  }

  comm_free(&dup);

  return 0;
}

#undef min
