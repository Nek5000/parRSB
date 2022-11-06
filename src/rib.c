#include "parrsb-impl.h"
#include "sort.h"

extern unsigned get_rsb_bin(uint id, uint np);

extern int power_serial(double *y, int N, double *A, int verbose);

static void calc_inertia(char *elems, uint nel, size_t usize, int ndim,
                         struct comm *c) {
  double avg[3] = {0.0, 0.0, 0.0};
  for (uint i = 0; i < nel; i++) {
    struct rcb_element *ei = (struct rcb_element *)(elems + i * usize);
    avg[0] += ei->coord[0], avg[1] += ei->coord[1], avg[2] += ei->coord[2];
  }

  double wrk[9];
  slong nelg = nel;
  if (c != NULL) {
    comm_allreduce(c, gs_double, gs_add, avg, 3, wrk);
    comm_allreduce(c, gs_long, gs_add, &nelg, 1, wrk);
  }

  avg[0] /= nelg, avg[1] /= nelg, avg[2] /= nelg;

  double I[3][3];
  for (uint i = 0; i < 3; i++)
    I[i][0] = I[i][1] = I[i][2] = 0.0;

  double x, y, z;
  for (uint i = 0; i < nel; i++) {
    struct rcb_element *ei = (struct rcb_element *)(elems + i * usize);
    x = ei->coord[0] - avg[0];
    y = ei->coord[1] - avg[1];
    z = ei->coord[2] - avg[2];
    I[0][0] += x * x, I[0][1] += x * y, I[0][2] += x * z;
    I[1][0] += y * x, I[1][1] += y * y, I[1][2] += y * z;
    I[2][0] += z * x, I[2][1] += z * y, I[2][2] += z * z;
  }

  if (c != NULL)
    comm_allreduce(c, gs_double, gs_add, I, 9, wrk);

  // FIXME: 2D does not work. ev[2] = 0 if 2D
  double ev[3];
  power_serial(ev, ndim, (double *)I, 0);

  for (uint i = 0; i < nel; i++) {
    struct rcb_element *ei = (struct rcb_element *)(elems + i * usize);
    x = ei->coord[0] - avg[0];
    y = ei->coord[1] - avg[1];
    z = ei->coord[2] - avg[2];
    ei->fiedler = x * ev[0] + y * ev[1] + z * ev[2];
  }
}

void rib_local(struct array *a, size_t usize, uint sidx, uint eidx, int ndim,
               buffer *bfr) {
  sint size = eidx - sidx;
  if (size <= 1)
    return;

  char *st = (char *)a->ptr + usize * sidx;
  calc_inertia(st, size, usize, ndim, NULL);

  if (usize == sizeof(struct rcb_element))
    sarray_sort(struct rcb_element, st, size, fiedler, 3, bfr);
  else if (usize == sizeof(struct rsb_element))
    sarray_sort(struct rsb_element, st, size, fiedler, 3, bfr);

  uint midx = (sidx + eidx) / 2;
  rib_local(a, usize, sidx, midx, ndim, bfr);
  rib_local(a, usize, midx, eidx, ndim, bfr);
}

int rib(struct array *elems, size_t usize, int ndim, struct comm *ci,
        buffer *bfr) {
  struct comm c;
  comm_dup(&c, ci);

  while (c.np > 1) {
    calc_inertia((char *)elems->ptr, elems->n, usize, ndim, &c);

    if (usize == sizeof(struct rcb_element))
      parallel_sort(struct rcb_element, elems, fiedler, gs_double, 0, 1, &c,
                    bfr);
    else if (usize == sizeof(struct rsb_element))
      parallel_sort(struct rsb_element, elems, fiedler, gs_double, 0, 1, &c,
                    bfr);

    unsigned bin = get_rsb_bin(c.id, c.np);

    struct comm t;
    comm_split(&c, bin, c.id, &t);
    comm_free(&c);
    comm_dup(&c, &t);
    comm_free(&t);
  }
  comm_free(&c);

  rib_local(elems, usize, 0, elems->n, ndim, bfr);

  return 0;
}
