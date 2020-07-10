#include "genmap-impl.h"

#include <math.h>
#include <stdio.h>

int GenmapGetInitVector(GenmapHandle h,GenmapComm c,int global,
    GenmapVector *initVec_)
{
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapCreateVector(initVec_,lelt);
  GenmapVector initVec=*initVec_;

  GenmapInt i;
  GenmapElements elements = GenmapGetElements(h);
  if(global > 0) {
#if defined(GENMAP_TQLI)
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = GenmapGetLocalStartIndex(h) + i + 1;
    }
#else
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = (GenmapScalar) elements[i].globalId;
    }
#endif
  } else {
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = elements[i].fiedler;
    }
  }

  GenmapOrthogonalizebyOneVector(h, c, initVec,
      GenmapGetNGlobalElements(h));
  GenmapScalar rtr = GenmapDotVector(initVec, initVec);
  GenmapGop(c, &rtr, 1, GENMAP_SCALAR, GENMAP_SUM);
  GenmapScalar rni = 1.0 / sqrt(rtr);
  GenmapScaleVector(initVec, initVec, rni);

  return 0;
}

void GenmapRSB(GenmapHandle h) {
  GenmapInt i;
  GenmapElements e = GenmapGetElements(h);
  GenmapScan(h, GenmapGetLocalComm(h));
  for(i = 0; i < GenmapGetNLocalElements(h); i++) {
    e[i].globalId  = GenmapGetLocalStartIndex(h) + i + 1;
    e[i].globalId0 = GenmapGetLocalStartIndex(h) + i + 1;
  }

  if(GenmapCommRank(GenmapGetGlobalComm(h)) == 0 && h->dbgLevel > 0)
    printf("running RSB ...\n"), fflush(stdout);

  crystal_init(&(h->cr), &(h->local->gsc));
  buffer buf0 = null_buffer;

  int level=0;

  while(GenmapCommSize(GenmapGetLocalComm(h)) > 1) {
#if defined(GENMAP_TQLI)
    int global = 1;
#else
    int global=(GenmapCommSize(GenmapGetLocalComm(h))== \
        GenmapCommSize(GenmapGetGlobalComm(h)));
#endif

    GenmapVector initVec;
    GenmapComm c=GenmapGetLocalComm(h);

    GenmapInt totalIter;

    int iter;
    int ipass = 0;
    int npass = 50;
    int maxIter = 50;
#if defined(GENMAP_LANCZOS)
    do {
      GenmapGetInitVector(h,c,global,&initVec);
      iter = GenmapFiedlerLanczos(h,c,initVec,maxIter,global);
      ipass++;
      global = 0;
      GenmapDestroyVector(initVec);
    } while(ipass < npass && iter == maxIter);

    totalIter=(ipass-1)*maxIter+iter;
#else
    maxIter=1000;
    GenmapGetInitVector(h,c,global,&initVec);
    totalIter=GenmapFiedlerPowerIter(h,c,initVec,maxIter,global);
    GenmapDestroyVector(initVec);
#endif

    GenmapInt iterMax,iterMin,iterSum;
    GenmapComm globalComm=GenmapGetGlobalComm(h);
    GenmapReduce(globalComm,&iterMax,&totalIter,1,GENMAP_INT,GENMAP_MAX);
    GenmapReduce(globalComm,&iterMin,&totalIter,1,GENMAP_INT,GENMAP_MIN);
    GenmapReduce(globalComm,&iterSum,&totalIter,1,GENMAP_INT,GENMAP_SUM);

    if(GenmapCommRank(GenmapGetGlobalComm(h))==0 && h->dbgLevel>1){
      printf("level %02d: %5d,%5d,%5d (min,mean,max) iterations\n",
          level,iterMax,iterSum/GenmapCommSize(globalComm),iterMin);
      fflush(stdout);
    }

    GenmapBinSort(h, GENMAP_FIEDLER, &buf0);

    int bin;
    GenmapInt np = GenmapCommSize(GenmapGetLocalComm(h));
    GenmapInt id = GenmapCommRank(GenmapGetLocalComm(h));
    if(id<(np+1)/2) bin = 0;
    else bin = 1;

    c = GenmapGetLocalComm(h);
    GenmapSplitComm(h, &c, bin);
    GenmapSetLocalComm(h, c);

#if defined(GENMAP_TQLI)
    GenmapBinSort(h, GENMAP_GLOBALID, &buf0);
#endif
    level++;
  }

  crystal_free(&(h->cr));
  buffer_free(&buf0);
}
