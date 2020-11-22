#include <genmap-impl.h>

void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

void matrix_inverse(int N,double *A){
  int size=N*N;
  int info;

  int *ipiv=(int*) calloc(N,sizeof(int));
  double *work=(double *) calloc(N*N,sizeof(double));

  dgetrf_(&N,&N,A,&N,ipiv,&info);
  if(info!=0){
    printf("dgetrf_: %d\n",info);
  }

  dgetri_(&N,A,&N,ipiv,work,&size,&info);
  if(info!=0){
    printf("dgetri_: %d\n",info);
  }

  free(ipiv);
  free(work);
}
