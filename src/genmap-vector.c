#include <genmap-impl.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

int genmap_vector_create(genmap_vector *x, GenmapInt size) {
  /* Asserts:
       - size > 0
  */
  assert(size > 0);

  GenmapMalloc(1, x);
  if (*x == NULL) {
    return 1;
  }

  (*x)->size = size;
  (*x)->data = NULL;

  GenmapMalloc((size_t)size, &(*x)->data);
  if ((*x)->data == NULL) {
    return 1;
  }

  return 0;
}

int genmap_destroy_vector(genmap_vector x) {
  if (x->data) {
    GenmapFree(x->data);
  }

  if (x) {
    GenmapFree(x);
  }

  return 0;
}

int genmap_vector_copy(genmap_vector y, genmap_vector x) {
  /* Asserts:
       - size y = size x
  */
  assert(y->size >= x->size);

  GenmapInt n = x->size;
  GenmapInt i;
  for (i = 0; i < n; i++) {
    y->data[i] = x->data[i];
  }

  return 0;
}

int GenmapCreateOnesVector(genmap_vector *x, GenmapInt size) {
  genmap_vector_create(x, size);

  GenmapInt i;
  for (i = 0; i < size; i++) {
    (*x)->data[i] = 1.;
  }

  return 0;
}

int genmap_vector_create_zeros(genmap_vector *x, GenmapInt size) {
  genmap_vector_create(x, size);

  GenmapInt i;
  for (i = 0; i < size; i++) {
    (*x)->data[i] = 0.;
  }

  return 0;
}
