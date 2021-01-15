#include <stdio.h>
#include <stdlib.h>

#include <genmap-impl.h>

int genmap_init(genmap_handle *h_, comm_ext ce, parRSB_options *options) {
  GenmapMalloc(1, h_);
  genmap_handle h = *h_;

  GenmapCreateComm(&h->global, ce);
  GenmapCreateComm(&h->local, ce);

  h->weights = NULL;

  h->options = options;

  return 0;
}

int genmap_finalize(genmap_handle h) {
  if (h->weights != NULL)
    GenmapFree(h->weights);
  if (GenmapGetGlobalComm(h))
    GenmapDestroyComm(GenmapGetGlobalComm(h));
  if (GenmapGetLocalComm(h))
    GenmapDestroyComm(h->local);

  GenmapFree(h);

  return 0;
}

int GenmapMallocArray(size_t n, size_t unit, void *p) {
  int ierr = posix_memalign((void **)p, GENMAP_ALIGN, n * unit);
  if (ierr)
    printf("GenmapMallocArray Failed: %s:%d\n", __FILE__, __LINE__);
  return ierr;
}

int GenmapCallocArray(size_t n, size_t unit, void *p) {
  int ierr = 0;
  *(void **)p = calloc(n, unit);
  if (n && unit && !*(void **)p) {
    ierr = 1;
    printf("GenmapCallocArray Failed: %s:%d\n", __FILE__, __LINE__);
  }
  return ierr;
}

int GenmapReallocArray(size_t n, size_t unit, void *p) {
  int ierr = 0;
  *(void **)p = realloc(*(void **)p, n * unit);
  if (n && unit && !*(void **)p) {
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
