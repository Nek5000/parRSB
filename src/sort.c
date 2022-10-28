#include "sort.h"
#include <float.h>
#include <math.h>
#include <stdarg.h>

struct sort {
  struct array *a;
  size_t usize, align;
  unsigned nfields;
  size_t offset[3];
  gs_dom t[3];
  buffer *buf;
};

//-----------------------------------------------------------------------------
// Common helper functions
//
static double get_scalar(struct array *a, uint i, uint offset, uint usize,
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

static void get_extrema(double extrema[2], struct sort *data, uint field,
                        struct comm *c) {
  struct array *a = data->a;
  uint usize = data->usize;
  uint offset = data->offset[field];
  gs_dom t = data->t[field];

  sint size = a->n;
  if (size == 0) {
    extrema[0] = -DBL_MAX;
    extrema[1] = -DBL_MAX;
  } else {
    extrema[0] = -get_scalar(a, 0, offset, usize, t);
    extrema[1] = get_scalar(a, size - 1, offset, usize, t);
  }

  double wrk[2];
  comm_allreduce(c, gs_double, gs_max, extrema, 2, wrk);
  extrema[0] *= -1;
}

static int set_dest(uint *proc, uint size, sint np, slong start, slong nelem) {
  if (nelem == 0)
    return 1;

  uint nelt = nelem / np, nrem = nelem - np * nelt;
  if (nrem == 0) {
    for (uint i = 0; i < size; i++) {
      proc[i] = (start + i) / nelt;
    }
  } else {
    uint s = np - nrem;
    slong t = nelt * s;
    for (uint i = 0; i < size; i++) {
      if (start + i < t)
        proc[i] = (start + i) / nelt;
      else
        proc[i] = s + (start + i - t) / (nelt + 1);
    }
  }

  return 0;
}

static int load_balance(struct array *a, size_t size, struct crystal *cr) {
  struct comm *c = &cr->comm;

  slong out[2][1], wrk[2][1], in = a->n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0], nelem = out[1][0];

  uint *proc = tcalloc(uint, a->n);
  set_dest(proc, a->n, c->np, start, nelem);
  sarray_transfer_ext_(a, size, proc, sizeof(uint), cr);
  free(proc);

  return 0;
}

//-----------------------------------------------------------------------------
// Parallel Bin-Sort
//
static int set_bin(uint **proc_, struct sort *s, uint field, struct comm *c) {
  struct array *a = s->a;
  gs_dom t = s->t[field];
  uint offset = s->offset[field];

  uint size = a->n;
  uint *proc = *proc_ = tcalloc(uint, size);

  double extrema[2];
  get_extrema(extrema, s, field, c);
  double range = extrema[1] - extrema[0];

  if (size == 0)
    return 0;

  sint np = c->np;
  uint id = 0;
  uint index = 0;
  do {
    double end = extrema[0] + (range / np) * (id + 1);
    while (index < size) {
      double val = get_scalar(a, index, offset, s->usize, t);
      if (val <= end)
        proc[index++] = id;
      else
        break;
    }
    id++;
  } while (id < np && index < size);
  for (; index < size; index++)
    proc[index] = np - 1;

  return 0;
}

static int sort_field(struct array *arr, size_t usize, gs_dom t, uint off,
                      buffer *buf, int keep) {
  uint nunits = arr->n;
  void *ptr = arr->ptr;
  switch (t) {
  case gs_double:
    gslib_sortp_double(buf, keep, (double *)((char *)ptr + off), nunits, usize);
    break;
  case gs_float:
    gslib_sortp_float(buf, keep, (float *)((char *)ptr + off), nunits, usize);
    break;
  case gs_long:
    gslib_sortp_ull(buf, keep, (ulong *)((char *)ptr + off), nunits, usize);
    break;
  case gs_int:
    gslib_sortp_ui(buf, keep, (uint *)((char *)ptr + off), nunits, usize);
    break;
  default:
    break;
  }

  return 0;
}

int sort_local(struct sort *s) {
  struct array *a = s->a;
  int i = s->nfields - 1;

  sort_field(a, s->usize, s->t[i], s->offset[i], s->buf, 0), i--;
  while (i >= 0)
    sort_field(a, s->usize, s->t[i], s->offset[i], s->buf, 1), i--;
  sarray_permute_buf_(s->align, s->usize, a->ptr, a->n, s->buf);

  return 0;
}

static int parallel_bin_sort(struct sort *s, struct comm *c) {
  // Local sort
  sort_local(s);

  // Set destination bin
  uint *proc;
  set_bin(&proc, s, 0, c);

  // Transfer to destination processor
  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer_ext_(s->a, s->usize, proc, sizeof(uint), &cr);
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

static int init_probes(struct hypercube *hcube, struct comm *c) {
  double extrema[2];
  get_extrema(extrema, hcube->data, 0, c);
  double range = extrema[1] - extrema[0];

  double delta = range / (hcube->nprobes - 1);
  hcube->probes[0] = extrema[0];
  hcube->probes[1] = extrema[0] + delta;
  hcube->probes[2] = extrema[1];

  return 0;
}

static int update_probe_counts(struct hypercube *hcube, struct comm *c) {
  for (uint i = 0; i < hcube->nprobes; i++)
    hcube->probe_cnt[i] = 0;

  struct sort *s = hcube->data;
  struct array *arr = s->a;
  for (uint e = 0; e < arr->n; e++) {
    double val_e = get_scalar(arr, e, s->offset[0], s->usize, s->t[0]);
    for (uint i = 0; i < hcube->nprobes; i++) {
      if (val_e < hcube->probes[i])
        hcube->probe_cnt[i]++;
    }
  }

  ulong buf[3];
  comm_allreduce(c, gs_long, gs_add, hcube->probe_cnt, hcube->nprobes, buf);

  return 0;
}

static int update_probes(slong nelem, double *probes, ulong *probe_cnt,
                         uint threshold) {
  slong expected = nelem / 2, cnt = probe_cnt[1];
  if (llabs(expected - cnt) < threshold)
    return 0;

  if (probe_cnt[1] < expected)
    probes[0] = probes[1];
  else
    probes[2] = probes[1];

  probes[1] = probes[0] + (probes[2] - probes[0]) / 2.0;

  return 0;
}

static int transfer_elem(struct hypercube *data, struct comm *c) {
  struct sort *input = data->data;
  uint usize = input->usize, offset = input->offset[0];
  gs_dom t = input->t[0];
  struct array *a = input->a;

  uint size = a->n, lown = 0, uppern = 0;
  for (uint e = 0; e < size; e++) {
    double val = get_scalar(a, e, offset, usize, t);
    if (val < data->probes[1])
      lown++;
    else
      uppern++;
  }

  slong out[2][2], in[2] = {lown, uppern}, buf[2][2];
  comm_scan(out, c, gs_long, gs_add, in, 2, buf);
  slong lstart = out[0][0], ustart = out[0][1];
  slong lelem = out[1][0], uelem = out[1][1];

  uint np = c->np, lnp = np / 2;
  uint *proc = tcalloc(uint, size);
  set_dest(proc, lnp, lstart, lown, lelem);
  set_dest(proc + lown, np - lnp, ustart, uppern, uelem);

  for (uint e = lown; e < size; e++)
    proc[e] += lnp;

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer_ext_(a, usize, proc, sizeof(uint), &cr);
  crystal_free(&cr);

  free(proc);

  return 0;
}

static int parallel_hcube_sort_recursive(struct hypercube *data,
                                         struct comm *c) {
  struct sort *input = data->data;
  struct array *a = input->a;
  gs_dom t = input->t[0];
  uint offset = input->offset[0];

  slong out[2][1], wrk[2][1], in = a->n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0], nelem = out[1][0];

  sint size = c->np, rank = c->id;
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
  slong cnt = data->probe_cnt[1];
  while (llabs(nelem / 2 - cnt) > threshold && iter++ < max_iter) {
    update_probes(nelem, data->probes, data->probe_cnt, threshold);
    update_probe_counts(data, c);
    cnt = data->probe_cnt[1];
  }

  transfer_elem(data, c);

  // Split the communicator and recurse
  struct comm nc;
  int bin = (rank < size / 2) ? 1 : 0;
  comm_split(c, bin, c->id, &nc);
  parallel_hcube_sort_recursive(data, &nc);
  comm_free(&nc);

  return 0;
}

static int parallel_hcube_sort(struct sort *s, struct comm *c) {
  struct hypercube hcube = {.data = s, .nprobes = 3};

  hcube.probes = tcalloc(double, hcube.nprobes);
  hcube.probe_cnt = tcalloc(ulong, hcube.nprobes);

  parallel_hcube_sort_recursive(&hcube, c);

  free(hcube.probes);
  free(hcube.probe_cnt);

  return 0;
}

int parallel_sort_(struct array *array, unsigned algo, unsigned balance,
                   size_t usize, size_t align, struct comm *ci, buffer *bfr,
                   unsigned nfields, ...) {
  struct sort s;

  s.a = array;
  s.usize = usize, s.align = align;
  s.nfields = nfields;

  va_list args;
  va_start(args, nfields);
  for (unsigned i = 0; i < nfields; i++) {
    s.offset[i] = va_arg(args, size_t);
    s.t[i] = va_arg(args, gs_dom);
  }
  va_end(args);

  s.buf = bfr;

  struct comm c;
  comm_dup(&c, ci);
  switch (algo) {
  case 0:
    parallel_bin_sort(&s, &c);
    break;
  case 1:
    parallel_hcube_sort(&s, &c);
    break;
  default:
    break;
  }
  comm_free(&c);

  if (balance) {
    struct crystal cr;
    crystal_init(&cr, ci);
    load_balance(array, usize, &cr);
    crystal_free(&cr);
    sort_local(&s);
  }

  return 0;
}
