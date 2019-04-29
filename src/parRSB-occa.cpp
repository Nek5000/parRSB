/*
 * Dedicated to all the victims of 4/21 - Easter day attack in SL.
*/

#include "parRSB-occa.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int parRSBOccaSetup(GenmapHandle h) {
  // OCCA build stuff
  char deviceConfig[BUFSIZ];

  GenmapComm c = GenmapGetLocalComm(h);
  int rank, size;
  size = GenmapCommSize(c);
  rank = GenmapCommRank(c);

  long int hostId = gethostid();
  long int* hostIds = (long int*) calloc(size, sizeof(long int));
  MPI_Allgather(&hostId, 1, MPI_LONG, hostIds, 1, MPI_LONG, GenmapGetMPIComm(c));

  int device_id = 0;
  int totalDevices = 0;
  for(int r = 0; r < rank; r++) {
    if(hostIds[r] == hostId) device_id++;
  }
  for(int r = 0; r < size; r++) {
    if(hostIds[r] == hostId) totalDevices++;
  }

  // read thread model/device/platform from options
  sprintf(deviceConfig, "mode: 'CUDA', device_id: %d", device_id);

  parRSBKrylov krylov = parRSBGetKrylov(h);
  // there is a memory leak here, we don't delete it anywhere
  krylov = new parRSBKrylov_private[1];
  krylov->device.setup((std::string)deviceConfig);

  // TODO: These should not be hardcoded
  occa::properties prop;
  prop["defines"].asObject();
  prop["includes"].asArray();
  prop["header"].asArray();
  prop["flags"].asObject();

  prop["defines/" "dlong"      ] = "long";
  prop["defines/" "dfloat"     ] = "double";
  prop["defines/" "p_blockSize"] = 256;

  krylov->innerProductKernel = krylov->device.buildKernel(
                                 PARRSB_OKL_DIR "/innerProduct.okl", "innerProduct", prop);
  krylov->weightedInnerProduct1Kernel = krylov->device.buildKernel(
                                          PARRSB_OKL_DIR "/weightedInnerProduct1.okl",
                                          "weightedInnerProduct1", prop);
  krylov->weightedInnerProduct2Kernel = krylov->device.buildKernel(
                                          PARRSB_OKL_DIR "/weightedInnerProduct2.okl",
                                          "weightedInnerProduct2", prop);
  krylov->scaledAddKernel = krylov->device.buildKernel(
                              PARRSB_OKL_DIR "/scaledAdd.okl", "scaledAdd", prop);
  krylov->dotMultiplyKernel = krylov->device.buildKernel(
                                PARRSB_OKL_DIR "/dotMultiply.okl", "dotMultiply", prop);
  krylov->dotMultiplyAddKernel = krylov->device.buildKernel(
                                   PARRSB_OKL_DIR "/dotMultiplyAdd.okl", "dotMultiplyAdd", prop);
  krylov->dotDivideKernel = krylov->device.buildKernel(PARRSB_OKL_DIR "/dotDivide.okl",
                              "dotDivide", prop);
  krylov->weightedNorm2Kernel = krylov->device.buildKernel(PARRSB_OKL_DIR "/weightedNorm2.okl",
                                  "weightedNorm2", prop);
  krylov->norm2Kernel = krylov->device.buildKernel(PARRSB_OKL_DIR "/norm2.okl", "norm2", prop);
}

int parRSBLaplacianSetup(GenmapHandle h) {
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapInt nv = GenmapGetNVertices(h);
  GenmapUInt numPoints = (GenmapUInt) nv * lelt;

  GenmapLong *vertices;
  parRSBGetVertices(h, vertices);

  parRSBKrylov krylov = parRSBGetKrylov(h);

  krylov->r = (GenmapScalar *) calloc(lelt, sizeof(GenmapScalar));
  krylov->p = (GenmapScalar *) calloc(lelt, sizeof(GenmapScalar));
  krylov->w = (GenmapScalar *) calloc(lelt, sizeof(GenmapScalar));
  krylov->weights = (GenmapScalar *) calloc(lelt, sizeof(GenmapScalar));

  krylov->o_r       = krylov->device.malloc(lelt * sizeof(GenmapScalar));
  krylov->o_p       = krylov->device.malloc(lelt * sizeof(GenmapScalar));
  krylov->o_w       = krylov->device.malloc(lelt * sizeof(GenmapScalar));
  krylov->o_weights = krylov->device.malloc(lelt * sizeof(GenmapScalar));
}
