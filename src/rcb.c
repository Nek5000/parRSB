#include "parrsb-impl.h"
#include "sort.h"

extern unsigned get_proc_bin(uint id, uint np);

static unsigned get_axis(size_t unit_size, char *elems, uint nel, int ndim,
                         struct comm *c) {
  double min[3] = {DBL_MAX, DBL_MAX, DBL_MAX},
         max[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};

  struct rcb_element *ei;
  for (uint i = 0; i < nel; i++) {
    ei = (struct rcb_element *)(elems + i * unit_size);
    if (ei->coord[0] < min[0])
      min[0] = ei->coord[0];
    if (ei->coord[0] > max[0])
      max[0] = ei->coord[0];

    if (ei->coord[1] < min[1])
      min[1] = ei->coord[1];
    if (ei->coord[1] > max[1])
      max[1] = ei->coord[1];
  }

  if (ndim == 3) {
    for (uint i = 0; i < nel; i++) {
      ei = (struct rcb_element *)(elems + i * unit_size);
      if (ei->coord[2] < min[2])
        min[2] = ei->coord[2];
      if (ei->coord[2] > max[2])
        max[2] = ei->coord[2];
    }
  }

  if (c != NULL) {
    double wrk[3];
    comm_allreduce(c, gs_double, gs_min, min, 3, wrk);
    comm_allreduce(c, gs_double, gs_max, max, 3, wrk);
  }

  unsigned axis = 0;
  for (unsigned d = 1; d < ndim; d++) {
    if ((max[axis] - min[axis]) < (max[d] - min[d]))
      axis = d;
  }

  return axis;
}

void rcb_local(struct array *a, size_t usize, uint sidx, uint eidx, int ndim,
               buffer *buf) {
  if (eidx <= sidx + 1)
    return;

  char *st = (char *)a->ptr + usize * sidx;

  uint size = eidx - sidx;
  unsigned axis = get_axis(usize, st, size, ndim, NULL);
  if (usize == sizeof(struct rcb_element)) {
    switch (axis) {
    case 0:
      sarray_sort(struct rcb_element, st, size, coord[0], 3, buf);
      break;
    case 1:
      sarray_sort(struct rcb_element, st, size, coord[1], 3, buf);
      break;
    case 2:
      sarray_sort(struct rcb_element, st, size, coord[2], 3, buf);
      break;
    default:
      break;
    }
  } else if (usize == sizeof(struct rsb_element)) {
    switch (axis) {
    case 0:
      sarray_sort(struct rsb_element, st, size, coord[0], 3, buf);
      break;
    case 1:
      sarray_sort(struct rsb_element, st, size, coord[1], 3, buf);
      break;
    case 2:
      sarray_sort(struct rsb_element, st, size, coord[2], 3, buf);
      break;
    default:
      break;
    }
  }

  uint midx = (sidx + eidx) / 2;
  rcb_local(a, usize, sidx, midx, ndim, buf);
  rcb_local(a, usize, midx, eidx, ndim, buf);
}

int rcb(struct array *elems, size_t usize, int ndim, struct comm *ci,
        buffer *bfr) {
  struct comm c;
  comm_dup(&c, ci);

  while (c.np > 1) {
    unsigned axis = get_axis(usize, (char *)elems->ptr, elems->n, ndim, &c);
    if (usize == sizeof(struct rcb_element)) {
      switch (axis) {
      case 0:
        parallel_sort(struct rcb_element, elems, coord[0], gs_double, 0, 1, &c,
                      bfr);
        break;
      case 1:
        parallel_sort(struct rcb_element, elems, coord[1], gs_double, 0, 1, &c,
                      bfr);
        break;
      case 2:
        parallel_sort(struct rcb_element, elems, coord[2], gs_double, 0, 1, &c,
                      bfr);
        break;
      default:
        break;
      }
    } else if (usize == sizeof(struct rsb_element)) {
      switch (axis) {
      case 0:
        parallel_sort(struct rsb_element, elems, coord[0], gs_double, 0, 1, &c,
                      bfr);
        break;
      case 1:
        parallel_sort(struct rsb_element, elems, coord[1], gs_double, 0, 1, &c,
                      bfr);
        break;
      case 2:
        parallel_sort(struct rsb_element, elems, coord[2], gs_double, 0, 1, &c,
                      bfr);
        break;
      default:
        break;
      }
    }

    unsigned bin = get_proc_bin(c.id, c.np);

    struct comm t;
    comm_split(&c, bin, c.id, &t);
    comm_free(&c);
    comm_dup(&c, &t);
    comm_free(&t);
  }
  comm_free(&c);

  rcb_local(elems, usize, 0, elems->n, ndim, bfr);

  return 0;
}
