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
  prop["defines/" "p_NV"] = 8;

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

  krylov->setScalar = krylov->device.buildKernel(
                              PARRSB_OKL_DIR "/setScalar.okl", "setScalar", prop);
  krylov->gatherFromVertices = krylov->device.buildKernel(
                              PARRSB_OKL_DIR "/setScalar.okl", "gatherFromVertices", prop);
  krylov->scatterToVertices = krylov->device.buildKernel(
                              PARRSB_OKL_DIR "/setScalar.okl", "scatterToVertices", prop);
  krylov->laplacian = krylov->device.buildKernel(
                              PARRSB_OKL_DIR "/setScalar.okl", "laplacian", prop);
}

int parRSBLaplacianSetup(GenmapHandle h, GenmapComm c) {
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapInt nv = GenmapGetNVertices(h);
  GenmapUInt numPoints = (GenmapUInt) nv * lelt;

  parRSBKrylov krylov = parRSBGetKrylov(h);

  GenmapLong *vertices;
  MPI_Comm comm = GenmapGetMPIComm(c);
  parRSBGetVertices(h, vertices);
  krylov->ogs = ogsSetup((dlong) numPoints, vertices, comm, 1, krylov->device);

  krylov->o_tmp = krylov->device.malloc(sizeof(GenmapScalar)*numPoints);
  krylov->o_weights = krylov->device.malloc(lelt * sizeof(GenmapScalar));

  krylov->setScalar(numPoints, 1.0, krylov->o_tmp);
  ogsGatherScatter(krylov->o_tmp, dfloatString, "add", krylov->ogs);
  krylov->gatherFromVertices(lelt, krylov->o_tmp, krylov->o_weights);

  //Todo free occa::memory
  GenmapFree(vertices);
}

int parRSBLaplacian(GenmapHandle h, GenmapComm c, occa::memory o_u, occa::memory o_v) {
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapInt nv = GenmapGetNVertices(h);
  GenmapUInt numPoints = (GenmapUInt) lelt*nv;

  parRSBKrylov krylov = parRSBGetKrylov(h);

  krylov->o_tmp = krylov->device.malloc(sizeof(GenmapScalar)*numPoints);
  krylov->scatterToVertices(lelt, o_u, krylov->o_tmp);

  ogsGatherScatter(krylov->o_tmp, dfloatString, "add", krylov->ogs);

  krylov->laplacian(lelt, o_u, krylov->o_weights, krylov->o_tmp, o_v);

  return 0;
}
