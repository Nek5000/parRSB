#include "genmap-impl.h"
#include "genmap-io.h"

#include <stdlib.h>
#include <stdio.h>

int GenmapInit(GenmapHandle *h, GenmapCommExternal ce) {
  GenmapMalloc(1, h);
  GenmapHandle h_ = *h;

  GenmapCreateComm(&h_->global, ce);
  GenmapCreateComm(&h_->local, ce);

  h_->elementArray.ptr = NULL;
  h_->elementArray.n = (*h)->elementArray.max = 0;

  h_->dbgLevel = 0;
  h_->printStat = 0;

  GenmapMalloc(1,&h_->histogram);
  return 0;
}

int GenmapFinalize(GenmapHandle h) {
  if(GenmapGetGlobalComm(h))
    GenmapDestroyComm(GenmapGetGlobalComm(h));
  if(GenmapGetLocalComm(h))
    GenmapDestroyComm(h->local);

  array_free(&(h->elementArray));

  GenmapFree(h->histogram);

  GenmapFree(h);

  return 0;
}

GenmapElements GenmapGetElements(GenmapHandle h) {
  return (GenmapElements) h->elementArray.ptr;
}

void GenmapSetElements(GenmapHandle h, GenmapElements elements) {
  h->elementArray.ptr = elements;
}

GenmapInt GenmapGetNLocalElements(GenmapHandle h) {
  return h->elementArray.n;
}

void GenmapSetNLocalElements(GenmapHandle h, GenmapInt localElements) {
  array_init(struct GenmapElement_private, &h->elementArray, localElements);
  h->elementArray.n = localElements;
}

GenmapLong GenmapGetNGlobalElements(GenmapHandle h) {
  return h->nel;
}

void GenmapSetNGlobalElements(GenmapHandle h, GenmapLong globalElements) {
  h->nel = globalElements;
}

GenmapLong GenmapGetLocalStartIndex(GenmapHandle h) {
  return h->start;
}

void GenmapSetLocalStartIndex(GenmapHandle h, GenmapLong localStart) {
  h->start = localStart;
}

int GenmapGetNVertices(GenmapHandle h) {
  return h->nv;
}

void GenmapSetNVertices(GenmapHandle h, int nVertices) {
  h->nv = nVertices;
}

void GenmapScan(GenmapHandle h, GenmapComm c) {
  GenmapLong out[2][1], buf[2][1];
  GenmapLong lelt = GenmapGetNLocalElements(h);
  comm_scan(out, &(c->gsComm), genmap_gs_long, gs_add, &lelt, 1, buf);
  GenmapSetLocalStartIndex(h, out[0][0]);
  GenmapSetNGlobalElements(h, out[1][0]);
}

int GenmapMallocArray(size_t n, size_t unit, void *p) {
  int ierr = posix_memalign((void **)p, GENMAP_ALIGN, n * unit);
  if(ierr)
    printf("GenmapMallocArray Failed: %s:%d\n", __FILE__, __LINE__);
  return ierr;
}

int GenmapCallocArray(size_t n, size_t unit, void *p) {
  int ierr = 0;
  *(void **)p = calloc(n, unit);
  if(n && unit && !*(void **)p) {
    ierr = 1;
    printf("GenmapCallocArray Failed: %s:%d\n", __FILE__, __LINE__);
  }
  return ierr;
}

int GenmapReallocArray(size_t n, size_t unit, void *p) {
  int ierr = 0;
  *(void **)p = realloc(*(void **)p, n * unit);
  if(n && unit && !*(void **)p) {
    ierr = 1;
    printf("GenmapReallocArray Failed: %s:%d\n", __FILE__, __LINE__);
  }
  return ierr;
}

int GenmapFree(void *p) {
  free(p);
  p = NULL;
  return 0;
}
