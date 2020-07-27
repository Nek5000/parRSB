#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

#include <genmap-impl.h>

void GenmapRSB(GenmapHandle h){
  int maxIter=50;
  int npass  =50;

  GenmapInt i;
  GenmapElements e = GenmapGetElements(h);
  GenmapScan(h, GenmapGetLocalComm(h));
  for(i = 0; i < GenmapGetNLocalElements(h); i++) {
    e[i].globalId =GenmapGetLocalStartIndex(h)+i+1;
    e[i].globalId0=GenmapGetLocalStartIndex(h)+i+1;
  }

  if(GenmapCommRank(GenmapGetGlobalComm(h)) == 0 && h->dbgLevel > 0)
    printf("running RSB "), fflush(stdout);

  crystal_init(&(h->cr), &(h->local->gsc));
  buffer buf0 = null_buffer;

  while(GenmapCommSize(GenmapGetLocalComm(h)) > 1) {
    if(GenmapCommRank(GenmapGetGlobalComm(h)) == 0
        && h->dbgLevel > 1) printf("."), fflush(stdout);

    GenmapComm local_c=GenmapGetLocalComm(h);
    GenmapInt np=GenmapCommSize(local_c);

#if defined(GENMAP_PAUL)
    int global=1;
#else
    int global=(np==GenmapCommSize(GenmapGetGlobalComm(h)));
#endif

    int ipass = 0;
    int iter;
    do {
      iter=GenmapFiedler(h,local_c,maxIter,global);
      ipass++;
      global = 0;
    } while(ipass < npass && iter == maxIter);

    GenmapBinSort(h, GENMAP_FIEDLER, &buf0);

    int bin;
    GenmapInt id=GenmapCommRank(local_c);
    if(id < (np + 1) / 2) bin = 0;
    else bin = 1;

    GenmapSplitComm(h,&local_c,bin);
    GenmapSetLocalComm(h,local_c);

#if defined(GENMAP_PAUL)
    GenmapBinSort(h,GENMAP_GLOBALID,&buf0);
#endif
  }

  crystal_free(&(h->cr));
  buffer_free(&buf0);
}
