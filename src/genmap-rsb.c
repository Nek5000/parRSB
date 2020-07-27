#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

#include <genmap-impl.h>

void GenmapRSB(GenmapHandle h,int verbose){
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

  int rank=GenmapCommRank(GenmapGetGlobalComm(h));
  int level=0;

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

    int ntot=0,ipass=0;
    int iter;
    do {
      //iter=GenmapFiedlerLanczos(h,local_c,maxIter,global);
      iter=GenmapFiedlerRQI(h,local_c,maxIter,global);
      ntot+=iter;
      global=0;
    }while(++ipass < npass && iter==maxIter);

    int min,max,sum,buf[1];
    min=ntot;
    comm_allreduce(&local_c->gsc,gs_int,gs_min,&min,1,buf); // min
    max=ntot;
    comm_allreduce(&local_c->gsc,gs_int,gs_max,&max,1,buf); // max
    sum=ntot;
    comm_allreduce(&local_c->gsc,gs_int,gs_add,&sum,1,buf); // sum

    if(rank==0 && verbose)
      printf("level=%02d, iter: min=%d max=%d avg=%d\n",
        level,min,max,sum/np);

    GenmapBinSort(h, GENMAP_FIEDLER, &buf0);

    GenmapInt id=GenmapCommRank(local_c);
    int bin;
    if(id<(np+1)/2) bin=0;
    else bin=1;

    GenmapSplitComm(h,&local_c,bin);
    GenmapSetLocalComm(h,local_c);

#if defined(GENMAP_PAUL)
    GenmapBinSort(h,GENMAP_GLOBALID,&buf0);
#endif

    level++;
  }

  crystal_free(&(h->cr));
  buffer_free(&buf0);
}
