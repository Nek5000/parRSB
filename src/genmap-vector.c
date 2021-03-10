#include <genmap-impl.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

int GenmapCreateVector(genmap_vector *x, GenmapInt size) {
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

int GenmapSetVector(genmap_vector x, GenmapScalar *array) {
  memcpy(x->data, array, sizeof(GenmapScalar) * (size_t)x->size);
  return 0;
}

int GenmapDestroyVector(genmap_vector x) {
  if (x->data) {
    GenmapFree(x->data);
  }

  if (x) {
    GenmapFree(x);
  }

  return 0;
}

int GenmapCopyVector(genmap_vector y, genmap_vector x) {
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

int GenmapScaleVector(genmap_vector y, genmap_vector x, GenmapScalar alpha) {
  /* asserts:
       - size x = size y
  */
  assert(x->size == y->size);

  GenmapInt n = x->size;
  GenmapInt i;
  for (i = 0; i < n; i++) {
    y->data[i] = alpha * x->data[i];
  }

  return 0;
}

int GenmapCreateOnesVector(genmap_vector *x, GenmapInt size) {
  GenmapCreateVector(x, size);

  GenmapInt i;
  for (i = 0; i < size; i++) {
    (*x)->data[i] = 1.;
  }

  return 0;
}

int GenmapCreateZerosVector(genmap_vector *x, GenmapInt size) {
  GenmapCreateVector(x, size);

  GenmapInt i;
  for (i = 0; i < size; i++) {
    (*x)->data[i] = 0.;
  }

  return 0;
}

GenmapScalar GenmapDotVector(genmap_vector y, genmap_vector x) {
  /* asserts:
       - size x = size y
  */
  assert(x->size == y->size);

  GenmapScalar result = 0.0;
  GenmapInt i;
  for (i = 0; i < x->size; i++) {
    result += x->data[i] * y->data[i];
  }

  return result;
}

int GenmapAxpbyVector(genmap_vector z, genmap_vector x, GenmapScalar alpha,
                      genmap_vector y, GenmapScalar beta) {
  assert(z->size == x->size);
  assert(z->size == y->size);

  GenmapInt n = z->size;
  GenmapInt i;
  for (i = 0; i < n; i++) {
    z->data[i] = alpha * x->data[i] + beta * y->data[i];
  }

  return 0;
}

int GenmapPrintVector(genmap_vector x) {
  /* Asserts:
       - size x > 0
  */
  assert(x->size > 0);

  printf("(%lf", x->data[0]);
  GenmapInt i;
  for (i = 1; i < x->size - 1; i++) {
    printf(", %.10lf", x->data[i]);
  }

  if (x->size > 1) {
    printf(", %.10lf)", x->data[x->size - 1]);
  } else {
    printf(")");
  }

  return 0;
}

/* Orthogonalize by 1-vector (vector of all 1's) */
int GenmapOrthogonalizebyOneVector(struct comm *c, genmap_vector q1,
                                   GenmapULong n) {
  GenmapInt i;
  GenmapScalar sum = 0.0, buf;
  for (i = 0; i < q1->size; i++)
    sum += q1->data[i];

  comm_allreduce(c, gs_double, gs_add, &sum, 1, &buf);
  sum /= n;

  for (i = 0; i < q1->size; i++)
    q1->data[i] -= sum;

  return 0;
}
