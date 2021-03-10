#include "genmap-impl.h"

struct rsb_element * genmap_get_elements(genmap_handle h) {
  return (struct rsb_element *)h->elements->ptr;
}
void genmap_set_elements(genmap_handle h, struct array *elements) {
  h->elements = elements;
}

int genmap_get_nvertices(genmap_handle h) { return h->nv; }
void genmap_set_nvertices(genmap_handle h, int nv) { h->nv = nv; }

GenmapInt genmap_get_nel(genmap_handle h) { return h->elements->n; }

GenmapULong genmap_get_partition_nel(genmap_handle h) { return h->nel; }
void genmap_set_partition_nel(genmap_handle h, GenmapULong globalElements) {
  h->nel = globalElements;
}

GenmapLong genmap_get_local_start_index(genmap_handle h) { return h->start; }
void genmap_set_local_start_index(genmap_handle h, GenmapLong localStart) {
  h->start = localStart;
}

void genmap_comm_scan(genmap_handle h, struct comm *c) {
  GenmapLong out[2][1], buf[2][1];
  GenmapLong lelt = genmap_get_nel(h);
  comm_scan(out, c, gs_long_long, gs_add, &lelt, 1, buf);
  genmap_set_local_start_index(h, out[0][0]);
  genmap_set_partition_nel(h, out[1][0]);
}
