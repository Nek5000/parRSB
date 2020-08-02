#ifndef _GENMAP_IMPL_H_
#define _GENMAP_IMPL_H_

#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#ifdef GENMAP_DEBUG
#include <stdio.h>
#endif

#include <genmap.h>
#include <genmap-multigrid-precon.h>

#define GENMAP_FIEDLER 0
#define GENMAP_GLOBALID 1
#define GENMAP_PROC 2
#define GENMAP_ORIGIN 3

struct GenmapComm_private {
  struct comm gsc;
  struct gs_data *gsh;
  csr_mat M;
  buffer buf;
  GenmapScalar *b;
};

struct GenmapElement_private {
  GenmapScalar fiedler;
  GenmapLong globalId;
  GenmapLong globalId0;
  GenmapLong vertices[8];
  GenmapInt proc;
  GenmapInt origin;
};

int GenmapCreateElements(GenmapElements *e);
int GenmapDestroyElements(GenmapElements e);
GenmapElements GenmapGetElements_default(GenmapHandle h);

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
  double time[16];
};

int GenmapCreateHandle(GenmapHandle h);
int GenmapDestroyHandle(GenmapHandle h);

struct GenmapVector_private {
  GenmapInt size;
  GenmapScalar *data;
};

#define GenmapMalloc(n, p) GenmapMallocArray ((n), sizeof(**(p)), p)
#define GenmapCalloc(n, p) GenmapCallocArray ((n), sizeof(**(p)), p)
#define GenmapRealloc(n, p) GenmapReallocArray((n), sizeof(**(p)), p)

void GenmapFiedlerMinMax(GenmapHandle h, GenmapScalar *min,
    GenmapScalar *max);
void GenmapGlobalIdMinMax(GenmapHandle h, GenmapLong *min,
    GenmapLong *max);
GenmapInt GenmapSetFiedlerBin(GenmapHandle h);
GenmapInt GenmapSetGlobalIdBin(GenmapHandle h);
void GenmapAssignBins(GenmapHandle h, int field, buffer *buf0);
void GenmapTransferToBins(GenmapHandle h, int field, buffer *buf0);
void GenmapBinSort(GenmapHandle h, int field, buffer *buf0);

#define MAXDIM 3
typedef struct{
  int proc;
  int orig;
  int seq;
  unsigned long long id;
  double coord[MAXDIM];
}elm_rcb;

int parRCB(struct comm *ci,struct array *a,int ndim);
void rcb_local(struct array *a,uint start,uint end,
    int ndim,buffer *buf);

#endif
