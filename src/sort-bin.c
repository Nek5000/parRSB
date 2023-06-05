#include "sort-impl.h"

static uint *set_proc_from_val(struct sort *s, uint field,
                               const struct comm *c) {
  struct array *a = s->a;
  gs_dom t = s->t[field];
  uint offset = s->offset[field];

  uint size = a->n;
  uint *proc = tcalloc(uint, size + 1);

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
        proc[index] = id, index++;
      else
        break;
    }
    id++;
  } while (id < np && index < size);
  for (; index < size; index++)
    proc[index] = np - 1;
  return 0;
}

void parallel_bin_sort(struct sort *s, const struct comm *c) {
  // Locally sort the array first.
  sort_local(s);

  // Set destination bin based on the field value.
  uint *proc = set_proc_from_val(s, 0, c);

  // Calculate the global array size. If it is zero, nothing to do, just return.
  struct array *arr = s->a;
  slong ng = arr->n, wrk[2];
  comm_allreduce(c, gs_long, gs_add, &ng, 1, wrk);
  if (ng == 0)
    return;

  // Initialize the crystal router.
  struct crystal cr;
  crystal_init(&cr, c);

  // Transfer the array elements to destination processor. To avoid message
  // sizes larger than INT_MAX, we calculate total message size and then figure
  // out how many transfers we need. Then we transfer array using that many
  // transfers.
  size_t usize = s->unit_size;
  uint nt = 2 * ((ng * usize + INT_MAX - 1) / INT_MAX);
  uint tsize = (arr->n + nt - 1) / nt;

  struct array brr, crr;
  array_init_(&brr, tsize + 1, usize, __FILE__, __LINE__);
  array_init_(&crr, arr->n + 1, usize, __FILE__, __LINE__);

  char *pe = (char *)arr->ptr;
  uint off = 0;
  for (unsigned i = 0; i < nt; i++) {
    // Copy from arr to brr.
    brr.n = 0;
    uint off1 = off + tsize;
    for (uint j = off; j < arr->n && j < off1; j++)
      array_cat_(usize, &brr, &pe[j * usize], 1, __FILE__, __LINE__);
    sarray_transfer_ext_(&brr, usize, &proc[off], sizeof(uint), &cr);
    array_cat_(usize, &crr, brr.ptr, brr.n, __FILE__, __LINE__);
    off = (off1 < arr->n ? off1 : arr->n);
  }
  array_free(&brr), free(proc);

  arr->n = 0;
  array_cat_(usize, arr, crr.ptr, crr.n, __FILE__, __LINE__);
  array_free(&crr);

  crystal_free(&cr);

  // Locally sort again to make sure that we have both globally and locally
  // sorted array.
  sort_local(s);
}
