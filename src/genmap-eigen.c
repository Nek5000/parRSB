#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <genmap-impl.h>

int GenmapSymTriDiagSolve(genmap_vector x, genmap_vector b, genmap_vector alpha,
                          genmap_vector beta) {
  assert((x->size == b->size) && (x->size == alpha->size));
  assert(alpha->size == beta->size + 1);
  assert(b->size > 0);

  GenmapInt n = b->size;

  genmap_vector diag;
  genmap_vector_create(&diag, n);
  genmap_vector_copy(diag, alpha);

  genmap_vector_copy(x, b);

  GenmapInt i;
  for (i = 0; i < n - 1; i++) {
    GenmapScalar m = (beta->data[i] / diag->data[i]);
    x->data[i + 1] = x->data[i + 1] - m * x->data[i];
    diag->data[i + 1] = diag->data[i + 1] - m * beta->data[i];
  }

  x->data[n - 1] = x->data[n - 1] / diag->data[n - 1];

  for (i = n - 2; i >= 0; i--) {
    x->data[i] = (x->data[i] - beta->data[i] * x->data[i + 1]) / diag->data[i];
  }

  genmap_destroy_vector(diag);
  return 0;
}

int genmap_power(double *y, int N, double *A, int verbose) {
  time_t t;
  srand((unsigned)time(&t));

  int i;
  GenmapScalar norm = 0.0;
  for (i = 0; i < N; i++) {
    y[i] = (rand() % 50) / 50.0;
    norm += y[i] * y[i];
  }

  GenmapScalar normi = 1.0 / sqrt(norm);
  for (i = 0; i < N; i++)
    y[i] *= normi;

  double *Ay;
  GenmapCalloc(N, &Ay);

  int j, k, l;
  GenmapScalar err = 1.0, lambda;
  for (i = 0; i < 100; i++) {
    norm = 0.0;
    for (j = 0; j < N; j++) {
      Ay[j] = 0.0;
      for (k = 0; k < N; k++) {
        Ay[j] += A[j * N + k] * y[k];
      }
      norm += Ay[j] * Ay[j];
    }

    if (i > 0)
      err = (sqrt(norm) - lambda) / lambda;
    lambda = sqrt(norm);

    normi = 1.0 / sqrt(norm);
    for (j = 0; j < N; j++)
      y[j] = Ay[j] * normi;

    if (fabs(err) < 1.e-14)
      break;
  }

  GenmapFree(Ay);

  return i;
}

#if defined(GENMAP_UNDERSCORE)
#define FNAME(x) TOKEN_PASTE(x, _)
#else
#define FNAME(x) x
#endif

#define FDGETRF FNAME(dgetrf)
#define FDGETRI FNAME(dgetri)

#if defined(GENMAP_BLAS)
void FDGETRF(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);
void FDGETRI(int *N, double *A, int *lda, int *IPIV, double *WORK, int *lwork,
             int *INFO);

void matrix_inverse(int N, double *A) {
  int size = N * N;
  int info;

  int *ipiv = (int *)calloc(N, sizeof(int));
  double *work = (double *)calloc(N * N, sizeof(double));

  FDGETRF(&N, &N, A, &N, ipiv, &info);
  if (info != 0)
    printf("dgetrf: %d\n", info);

  FDGETRI(&N, A, &N, ipiv, work, &size, &info);
  if (info != 0)
    printf("dgetri: %d\n", info);

  free(ipiv);
  free(work);
}
#else
void matrix_inverse(int N, double *A) {}
#endif // GENMAP_BLAS

#undef FDGETRF
#undef FDGETRI
#undef FNAME

int genmap_inverse_power(double *y, int N, double *A, int verbose) {
  double *Ainv;
  GenmapCalloc(N * N, &Ainv);

  int j, k;
  for (j = 0; j < N; j++) {
    for (k = 0; k < N; k++)
      Ainv[j * N + k] = A[k * N + j];
  }

  matrix_inverse(N, Ainv);

  for (j = 0; j < N; j++) {
    for (k = 0; k < N; k++)
      A[j * N + k] = Ainv[k * N + j];
  }

  j = genmap_power(y, N, Ainv, verbose);

  GenmapFree(Ainv);

  return j;
}
