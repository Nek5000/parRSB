#ifndef _GENMAP_IMPL_H_
#define _GENMAP_IMPL_H_

#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
//
// Header for GPU occa implementation
//
#if defined(PARRSB_GPU)
#  include <occa.h>
#  include <parRSB-occa.h>
#endif
//
// Header for gslib and genmap API
//
#include "genmap-gslib.h"
#include "genmap.h"

#ifdef __cplusplus
extern "C" {
#endif
//
// Fiedler fields
//
#define GENMAP_FIEDLER 0
#define GENMAP_GLOBALID 1
#define GENMAP_PROC 2
#define GENMAP_ORIGIN 3
//
// GenmapComm
//
struct GenmapComm_private {
  struct comm gsComm;
  struct gs_data *verticesHandle;
  GenmapScalar *laplacianWeights;
  buffer buf;
};
//
// GenmapElements
//
struct GenmapElement_private {
  GenmapScalar fiedler;
  GenmapLong globalId;
  GenmapLong globalId0;
  GenmapLong vertices[8];
  GenmapInt proc;
  GenmapInt origin;
};
//
// GenmapElements: Create, Destroy
//
int GenmapCreateElements(GenmapElements *e);
int GenmapDestroyElements(GenmapElements e);
GenmapElements GenmapGetElements_default(GenmapHandle h);
//
// GenmapHandle
//
struct GenmapHandle_private {
  GenmapComm global;
  GenmapComm local;

  GenmapLong nel;
  GenmapLong Nnodes;
  GenmapLong start;
  int nv;

  struct array elementArray;

  struct crystal cr;

  int dbgLevel;
  int printStat;

#if defined(PARRSB_GPU)
  parRSBKrylov krylov;
#endif
};
//
// GenmapHandle
//
int GenmapCreateHandle(GenmapHandle h);
int GenmapDestroyHandle(GenmapHandle h);
//
// GenmapVector
//
struct GenmapVector_private {
  GenmapInt size;
  GenmapScalar *data;
};

#define GenmapMalloc(n, p) GenmapMallocArray ((n), sizeof(**(p)), p)
#define GenmapCalloc(n, p) GenmapCallocArray ((n), sizeof(**(p)), p)
#define GenmapRealloc(n, p) GenmapReallocArray((n), sizeof(**(p)), p)

#endif

#ifdef __cplusplus
}
#endif
