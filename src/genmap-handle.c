#include "genmap-impl.h"

GenmapElements GenmapGetElements(genmap_handle h) {
  return (GenmapElements)h->elements->ptr;
}

GenmapInt GenmapGetNLocalElements(genmap_handle h) { return h->elements->n; }

void genmap_set_elements(genmap_handle h, struct array *localElements) {
  h->elements = localElements;
}

GenmapLong genmap_get_global_nel(genmap_handle h) { return h->nel; }

void GenmapSetNGlobalElements(genmap_handle h, GenmapLong globalElements) {
  h->nel = globalElements;
}

GenmapLong GenmapGetLocalStartIndex(genmap_handle h) { return h->start; }

void GenmapSetLocalStartIndex(genmap_handle h, GenmapLong localStart) {
  h->start = localStart;
}

int GenmapGetNVertices(genmap_handle h) { return h->nv; }

void genmap_set_vertices(genmap_handle h, int nVertices) { h->nv = nVertices; }

void genmap_scan(genmap_handle h, genmap_comm c) {
  GenmapLong out[2][1], buf[2][1];
  GenmapLong lelt = GenmapGetNLocalElements(h);
  comm_scan(out, &(c->gsc), gs_long_long, gs_add, &lelt, 1, buf);
  GenmapSetLocalStartIndex(h, out[0][0]);
  GenmapSetNGlobalElements(h, out[1][0]);
}

void genmap_comm_scan(genmap_handle h, struct comm *c) {
  GenmapLong out[2][1], buf[2][1];
  GenmapLong lelt = GenmapGetNLocalElements(h);
  comm_scan(out, c, gs_long_long, gs_add, &lelt, 1, buf);
  GenmapSetLocalStartIndex(h, out[0][0]);
  GenmapSetNGlobalElements(h, out[1][0]);
}
