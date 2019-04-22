/*
 * Dedicated to all the victims of 4/21 - Easter day attack in SL.
*/

#include "genmap-occa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int GenmapOccaDeviceConfig(GenmapHandle h) {
  // OCCA build stuff
  char deviceConfig[BUFSIZ];

  GenmapComm c = GenmapGetLocalComm(h);
  int rank, size;
  size = GenmapCommSize(c);
  rank = GenmapCommRank(c);

  long int hostId = gethostid();
  long int* hostIds = (long int*) calloc(size, sizeof(long int));
  MPI_Allgather(&hostId, 1, MPI_LONG, hostIds, 1, MPI_LONG, c->gsComm.c);

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

  if(h->dbgLevel > 0)
    printf("Rank %d: device_id = %d \n", rank, device_id);

  h->device = occaCreateDeviceFromString(deviceConfig);

  GenmapKrylov k; GenmapMalloc(1, &k);
  GenmapOccaLaplacianOperatorSetup(h, k);
}

int GenmapOccaLaplacianOperatorSetup(GenmapHandle h, GenmapKrylov krylov) {
  occaProperties prop = occaCreateProperties();

  // TODO: These should not be hardcoded
  occaPropertiesSet(prop, "defines/dlong", occaString("long"));
  occaPropertiesSet(prop, "defines/dfloat", occaString("double"));
  occaPropertiesSet(prop, "defines/p_blockSize", occaInt(256));

  GenmapInt nLocal = GenmapGetNLocalElements(h);
  krylov->r = (GenmapScalar *) calloc(nLocal, sizeof(GenmapScalar));
  krylov->p = (GenmapScalar *) calloc(nLocal, sizeof(GenmapScalar));
  krylov->w = (GenmapScalar *) calloc(nLocal, sizeof(GenmapScalar));
  krylov->weights = (GenmapScalar *) calloc(nLocal, sizeof(GenmapScalar));

  krylov->o_r = occaDeviceMalloc(h->device, nLocal * sizeof(GenmapScalar), NULL,
                                 occaDefault);
  krylov->o_p = occaDeviceMalloc(h->device, nLocal * sizeof(GenmapScalar), NULL,
                                 occaDefault);
  krylov->o_w = occaDeviceMalloc(h->device, nLocal * sizeof(GenmapScalar), NULL,
                                 occaDefault);
  krylov->o_weights = occaDeviceMalloc(h->device, nLocal * sizeof(GenmapScalar),
                                       NULL, occaDefault);

  krylov->innerProductKernel = occaDeviceBuildKernel(h->device,
                               PARRSB_OKL_DIR "/innerProduct.okl", "innerProduct", prop);
}
