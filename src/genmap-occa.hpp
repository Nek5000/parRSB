/*
 * Dedicated to all the victims of 4/21 - Easter day attack in SL.
*/

#ifndef _GENMAP_OCCA_HPP_
#define _GENMAP_OCCA_HPP_

#include "occa.hpp"

extern "C" {
  #include "genmap-occa.h"
  #include "genmap-impl.h"
}

struct GenmapKrylov_ {
  ogs_t *ogs;

  occa::memory o_r;
  occa::memory o_p;
  occa::memory o_w;
  occa::memory o_weights;

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
}; 

#endif
