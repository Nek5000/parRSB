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

  GenmapVector weights;
  struct array elementArray;
  struct crystal cr;

  int dbgLevel;
  int printStat;
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

/* Fiedler related */
void GenmapFiedlerMinMax(GenmapHandle h,GenmapScalar *min,
  GenmapScalar *max);
void GenmapGlobalIdMinMax(GenmapHandle h,GenmapLong *min,
  GenmapLong *max);
GenmapInt GenmapSetFiedlerBin(GenmapHandle h);
GenmapInt GenmapSetGlobalIdBin(GenmapHandle h);
void GenmapAssignBins(GenmapHandle h, int field, buffer *buf0);
void GenmapTransferToBins(GenmapHandle h, int field, buffer *buf0);
void GenmapBinSort(GenmapHandle h, int field, buffer *buf0);

/* Genmap Metrics */
typedef enum{
  ALLREDUCE,
  AXISLEN,
  BINN1,
  BINN2,
  COMMSPLIT,
  CRYSTAL,
  END,
  FIEDLER,
  FMG,
  GSOP,
  GSSETUP,
  LANCZOS,
  LAPLACIANSETUP,
  LAPLACIAN,
  LOCALSORT,
  LOADBALANCE0,
  LOADBALANCE1,
  NCONN,
  NNEIGHBORS,
  NFIEDLER,
  NDISCON,
  NPROJECTPF,
  PAIRWISE,
  PARSORT,
  PRECONSETUP,
  PRECONVCYCLE,
  PRECONAX,
  PROJECTPF,
  RCB,
  RCBTRANSFER,
  RSB,
  RQI,
  SETPROC,
  UPDATEPROBE,
} metric;

void metric_init();
void metric_finalize();
void metric_acc(metric m,double count);
void metric_tic(struct comm *c,metric m);
void metric_toc(struct comm *c,metric m);
double metric_get_value(int level,metric m);
void metric_push_level();
uint metric_get_levels();
void metric_print(struct comm *c);

/* genCon */
typedef struct{
  GenmapULong sequenceId;
  int nNeighbors;
  GenmapULong elementId;
  GenmapULong vertexId;
  uint workProc;
} vertex;

/* Components */
sint is_disconnected(struct comm *c,struct gs_data *gsh,buffer *buf,
  uint nelt,uint nv);

/* parRCB internals */
#define MAXDIM 3
typedef struct{
  int proc;
  int orig;
  int seq;
  unsigned long long id;
  double coord[MAXDIM];
} elm_rcb;

int parRCB(struct comm *ci,struct array *a,int ndim);
void rcb_local(struct array *a,uint start,uint end,int ndim,buffer *buf);

#endif
