#include <stdio.h>
#include <stdlib.h>

#include <genmap-impl.h>

int genmap_init(genmap_handle *h_, comm_ext ce, parrsb_options *options) {
  GenmapMalloc(1, h_);
  genmap_handle h = *h_;

  GenmapMalloc(1, &h->global);
  comm_init(h->global, ce);

  /* Weighted Laplacian */
  h->weights = NULL;
  h->gsw = NULL;
  buffer_init(&h->buf, 1024);

  h->options = options;

  return 0;
}

int genmap_finalize(genmap_handle h) {
  /* Weighted Laplacian */
  if (h->weights != NULL)
    GenmapFree(h->weights);
  if (h->gsw != NULL)
    gs_free(h->gsw);
  buffer_free(&h->buf);

  if (h->global != NULL) {
    comm_free(h->global);
    GenmapFree(h->global);
  }

  GenmapFree(h);

  return 0;
}

void *genmap_get_elements(genmap_handle h) {
  return (struct rsb_element *)h->elements->ptr;
}
int genmap_get_nvertices(genmap_handle h) { return h->nv; }
GenmapInt genmap_get_nel(genmap_handle h) { return h->elements->n; }
GenmapULong genmap_get_partition_nel(genmap_handle h) { return h->nel; }

void genmap_barrier(struct comm *c) {
#if defined(GENMAP_SYNC_BY_REDUCTION)
  sint dummy = c->id;
  sint buf;
  comm_allreduce(c, gs_int, gs_max, &dummy, 1, &buf);
#else
  comm_barrier(c);
#endif
}

void genmap_comm_split(struct comm *old, int bin, int key, struct comm *new_) {
  MPI_Comm new_comm;
  MPI_Comm_split(old->c, bin, key, &new_comm);
  comm_init(new_, new_comm);
  MPI_Comm_free(&new_comm);
}

void genmap_comm_scan(genmap_handle h, struct comm *c) {
  GenmapLong out[2][1], buf[2][1];
  GenmapLong lelt = genmap_get_nel(h);
  comm_scan(out, c, gs_long_long, gs_add, &lelt, 1, buf);
  h->nel = out[1][0];
  h->start = out[0][0];
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
