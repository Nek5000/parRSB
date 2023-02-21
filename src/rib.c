#include "parrsb-impl.h"
#include "sort.h"

extern int power_serial(double *y, int N, double *A, int verbose);
extern int inv_power_serial(double *y, uint N, double *A, int verbose);

static void get_rib_proj(char *elems, uint nel, size_t unit_size, int ndim,
                         struct comm *c) {
  double avg[3] = {0, 0, 0};
  for (uint i = 0; i < nel; i++) {
    struct rcb_element *ei = (struct rcb_element *)(elems + i * unit_size);
    avg[0] += ei->coord[0], avg[1] += ei->coord[1], avg[2] += ei->coord[2];
  }
  slong nelg = nel;

  double buf[9];
  if (c != NULL) {
    comm_allreduce(c, gs_double, gs_add, avg, 3, buf);
    comm_allreduce(c, gs_long, gs_add, &nelg, 1, buf);
  }

  avg[0] /= nelg, avg[1] /= nelg, avg[2] /= nelg;

  double I[9];
  for (unsigned i = 0; i < 9; i++)
    I[i] = 0;

  double x, y, z;
  for (uint i = 0; i < nel; i++) {
    struct rcb_element *ei = (struct rcb_element *)(elems + i * unit_size);
    x = ei->coord[0] - avg[0];
    y = ei->coord[1] - avg[1];
    z = ei->coord[2] - avg[2];
    I[0] += x * x, I[1] += x * y, I[2] += x * z;
    I[3] += y * x, I[4] += y * y, I[5] += y * z;
    I[6] += z * x, I[7] += z * y, I[8] += z * z;
  }

  if (c != NULL)
    comm_allreduce(c, gs_double, gs_add, I, 9, buf);

  // FIXME: 2D does not work
  double ev[3];
  power_serial(ev, ndim, (double *)I, 0);

  double norm = 0;
  for (unsigned i = 0; i < ndim; i++)
    norm += ev[i] * ev[i];
  norm = sqrt(norm);

  for (unsigned i = 0; i < ndim; i++)
    ev[i] /= norm;

  for (uint i = 0; i < nel; i++) {
    struct rcb_element *ei = (struct rcb_element *)(elems + i * unit_size);
    x = ei->coord[0] - avg[0];
    y = ei->coord[1] - avg[1];
    z = ei->coord[2] - avg[2];
    ei->fiedler = x * ev[0] + y * ev[1] + z * ev[2];
  }
}

void rib_local(struct array *a, size_t unit_size, uint start, uint end,
               int ndim, buffer *buf) {
  sint size = end - start;
  if (size <= 1)
    return;

  char *st = (char *)a->ptr + unit_size * start;
  get_rib_proj(st, size, unit_size, ndim, NULL);

  if (unit_size == sizeof(struct rcb_element))
    sarray_sort(struct rcb_element, st, size, fiedler, 3, buf);
  else if (unit_size == sizeof(struct rsb_element))
    sarray_sort(struct rsb_element, st, size, fiedler, 3, buf);

  uint mid = (start + end) / 2;
  rib_local(a, unit_size, start, mid, ndim, buf);
  rib_local(a, unit_size, mid, end, ndim, buf);
}

static int rib_level(struct array *a, size_t unit_size, int ndim,
                     struct comm *c, buffer *bfr) {
  if (c->np == 1)
    return 0;

  get_rib_proj((char *)a->ptr, a->n, unit_size, ndim, c);

  if (unit_size == sizeof(struct rcb_element))
    parallel_sort(struct rcb_element, a, fiedler, gs_double, 0, 1, c, bfr);
  else if (unit_size == sizeof(struct rsb_element))
    parallel_sort(struct rsb_element, a, fiedler, gs_double, 0, 1, c, bfr);

  return 0;
}

int rib(struct array *elements, size_t unit_size, int ndim, struct comm *ci,
        buffer *bfr) {
  struct comm c, t;
  comm_dup(&c, ci);

  int size = c.np, rank = c.id;
  while (size > 1) {
    rib_level(elements, unit_size, ndim, &c, bfr);

    int bin = 1;
    if (rank < (size + 1) / 2)
      bin = 0;

    comm_split(&c, bin, rank, &t);
    comm_free(&c);
    comm_dup(&c, &t);
    comm_free(&t);

    size = c.np, rank = c.id;
  }
  comm_free(&c);

  rib_local(elements, unit_size, 0, elements->n, ndim, bfr);

  return 0;
}
