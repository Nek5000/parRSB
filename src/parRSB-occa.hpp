#ifndef _PARRSB_OCCA_HPP_
#define _PARRSB_OCCA_HPP_

#include "genmap-types.h"

#include "occa.hpp"
#include "ogs.hpp"
#include "parRSB-occa.h"

struct parRSBKrylov_private {
  //
  // Device
  //
  occa::device device;
  //
  // Device and host memory pointers
  //
  occa::memory o_r;
  GenmapScalar *r;

  occa::memory o_p;
  GenmapScalar *p;

  occa::memory o_w;
  GenmapScalar *w;

  occa::memory o_weights;
  GenmapScalar *weights;
  //
  // Kernels
  //
  occa::kernel AxKernel;
  occa::kernel innerProductKernel;
  occa::kernel weightedInnerProduct1Kernel;
  occa::kernel weightedInnerProduct2Kernel;
  occa::kernel scaledAddKernel;
  occa::kernel dotMultiplyKernel;
  occa::kernel dotMultiplyAddKernel;
  occa::kernel dotDivideKernel;
  occa::kernel weightedNorm2Kernel;
  occa::kernel norm2Kernel;
  //
  // TODO: Add ogs handles
  //
  ogs_t *ogs;
};

#endif
