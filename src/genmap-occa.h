#ifdef __cplusplus
extern "C" {
#endif

#ifndef _GENMAP_OCCA_H_
#define _GENMAP_OCCA_H_

#include "occa.h"
#include "genmap-impl.h"

struct GenmapKrylov_private {
  //
  // Device and host memory pointers
  //
  occaMemory o_r;
  GenmapScalar *r;

  occaMemory o_p;
  GenmapScalar *p;

  occaMemory o_w;
  GenmapScalar *w;

  occaMemory o_weights;
  GenmapScalar *weights;
  //
  // Kernels
  //
  occaKernel AxKernel;
  occaKernel innerProductKernel;
  occaKernel weightedInnerProduct1Kernel;
  occaKernel weightedInnerProduct2Kernel;
  occaKernel scaledAddKernel;
  occaKernel dotMultiplyKernel;
  occaKernel dotMultiplyAddKernel;
  occaKernel dotDivideKernel;
  occaKernel weightedNorm2Kernel;
  occaKernel norm2Kernel;
  //
  // TODO: Add ogs handles
};

typedef struct GenmapKrylov_private* GenmapKrylov;

int GenmapOccaDeviceConfig(GenmapHandle h);
int GenmapOccaLaplacianOperatorSetup(GenmapHandle h, GenmapKrylov k);
int GenmapOccaLanczosSetup(GenmapHandle h, GenmapKrylov k);
int GenmapOccaLanczos(GenmapHandle h, GenmapKrylov k);

#endif

#ifdef __cplusplus
}
#endif
