#include "sort-impl.h"
#include <float.h>
#include <math.h>

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
                 const struct comm *c) {
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

  double buf[4];
  comm_allreduce(c, gs_double, gs_max, extrema, 2, buf);
  extrema[0] *= -1;
}

void set_proc_from_idx(uint *proc, uint size, sint np, slong start,
                       slong nelem) {
  if (nelem == 0)
    return;

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
  case gs_long: // FIXME gs_ulong
    gslib_sortp_ull(buf, keep, (ulong *)((char *)ptr + off), nunits, usize);
    break;
  case gs_int: // FIXME gs_uint
    gslib_sortp_ui(buf, keep, (uint *)((char *)ptr + off), nunits, usize);
    break;
  default:
    break;
  }

  return 0;
}

void sort_local(struct sort *s) {
  struct array *a = s->a;
  buffer *buf = s->buf;
  size_t usize = s->unit_size;
  int i = s->nfields - 1;

  sort_field(a, usize, s->t[i], s->offset[i], buf, 0), i--;
  while (i >= 0)
    sort_field(a, usize, s->t[i], s->offset[i], buf, 1), i--;
  sarray_permute_buf_(s->align, usize, a->ptr, a->n, buf);
}

static int load_balance(struct array *a, size_t size, const struct comm *c) {
  slong out[2][1], buf[2][1], in = a->n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
  slong start = out[0][0], nelem = out[1][0];

  uint *proc = tcalloc(uint, a->n + 1);
  set_proc_from_idx(proc, a->n, c->np, start, nelem);
  sarray_transfer_chunk(a, size, proc, c);
  free(proc);

  return 0;
}

void sarray_transfer_chunk(struct array *arr, const size_t usize,
                           const uint *proci, const struct comm *c) {
  // Calculate the global array size. If it is zero, nothing to do, just return.
  slong ng = arr->n, wrk[2];
  comm_allreduce(c, gs_long, gs_add, &ng, 1, wrk);
  if (ng == 0)
    return;

  // Initialize the crystal router.
  struct crystal cr;
  crystal_init(&cr, c);

  // Allocate `proc` with some buffer space.
  uint *proc = tcalloc(uint, arr->n + 1);
  for (uint i = 0; i < arr->n; i++)
    proc[i] = proci[i];

  // Transfer the array elements to destination processor. To avoid message
  // sizes larger than INT_MAX, we calculate total message size and then figure
  // out how many transfers we need. Then we transfer array using that many
  // transfers.
  slong msg_size = INT_MAX;
  uint nt = (ng * usize + msg_size - 1) / msg_size;
  uint tsize = (arr->n + nt - 1) / nt;

  struct array brr, crr;
  array_init_(&brr, tsize + 1, usize, __FILE__, __LINE__);
  array_init_(&crr, arr->n + 1, usize, __FILE__, __LINE__);

  char *pe = (char *)arr->ptr;
  uint off = 0;
  for (unsigned i = 0; i < nt; i++) {
    // Copy a chunk from `arr` to `brr`.
    brr.n = 0;
    uint off1 = off + tsize;
    for (uint j = off; j < arr->n && j < off1; j++)
      array_cat_(usize, &brr, &pe[j * usize], 1, __FILE__, __LINE__);
    assert(off <= arr->n);
    // Transfer the chunk in `brr` to the destination.
    sarray_transfer_ext_(&brr, usize, &proc[off], sizeof(uint), &cr);
    // Append the received chunk to `crr`.
    array_cat_(usize, &crr, brr.ptr, brr.n, __FILE__, __LINE__);
    off = (off1 < arr->n ? off1 : arr->n);
  }
  array_free(&brr);

  arr->n = 0;
  array_cat_(usize, arr, crr.ptr, crr.n, __FILE__, __LINE__);
  array_free(&crr);

  free(proc);
  crystal_free(&cr);
}

void parallel_sort_private(struct sort *data, const struct comm *c) {
  struct comm dup;
  comm_dup(&dup, c);

  int balance = data->balance, algo = data->algo;
  struct array *a = data->a;
  size_t usize = data->unit_size;

  struct hypercube hdata;
  switch (algo) {
  case 0:
    parallel_bin_sort(data, c);
    break;
  case 1:
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
    load_balance(a, usize, c);
    sort_local(data);
  }

  comm_free(&dup);
}
