#include "genmap-impl.h"

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

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

/* find the fiedler vector using Lanczos */
int GenmapFiedlerLanczos(GenmapHandle h, GenmapComm c,
    GenmapVector initVec,int maxIter, int global)
{
  GenmapVector alphaVec, betaVec;
  GenmapCreateVector(&alphaVec, maxIter);
  GenmapCreateVector(&betaVec, maxIter - 1);
  GenmapVector *q = NULL;

#if defined(GENMAP_TQLI)
  int iter = GenmapLanczosLegendary(h, c, initVec, maxIter, &q, alphaVec,
                                    betaVec);
#else
  int iter = GenmapLanczos         (h, c, initVec, maxIter, &q, alphaVec,
                                    betaVec);
#endif
  GenmapVector evLanczos, evTriDiag;
  GenmapCreateVector(&evTriDiag, iter);

  GenmapInt i;
#if defined(GENMAP_TQLI)
  /* Use TQLI and find the minimum eigenvalue and associated vector */
  GenmapVector *eVectors, eValues;
  GenmapTQLI(h, alphaVec, betaVec, &eVectors, &eValues);

  GenmapScalar eValMin = fabs(eValues->data[0]);
  GenmapInt eValMinI = 0;
  for(i = 1; i < iter; i++) {
    if(fabs(eValues->data[i]) < eValMin) {
      eValMin = fabs(eValues->data[i]);
      eValMinI = i;
    }
  }
  GenmapCopyVector(evTriDiag, eVectors[eValMinI]);
#else
  GenmapVector init;
  GenmapCreateVector(&init, iter);
  for(i = 0; i < iter; i++) {
    init->data[i] = i + 1.0;
  }
  GenmapScalar avg = 0.5 * iter * (1.0 + iter) / iter;
  for(i = 0; i < iter; i++) {
    init->data[i] -= avg;
  }
  GenmapInvPowerIter(evTriDiag, alphaVec, betaVec, init, 100);
  GenmapDestroyVector(init);
#endif

  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapCreateZerosVector(&evLanczos, lelt);
  GenmapInt j;
  for(i = 0; i < lelt; i++) {
    for(j = 0; j < iter; j++) {
      evLanczos->data[i] += q[j]->data[i] * evTriDiag->data[j];
    }
  }
#if defined(GENMAP_DEBUG)
  if(GenmapCommRank(GenmapGetGlobalComm(h)) == 0) {
    for(i = 0; i < evLanczos->size; i++) {
      printf("evLanczos:"GenmapScalarFormat"\n", evLanczos->data[i]);
    }
  }
#endif

  GenmapScalar lNorm = 0;
  for(i = 0; i < lelt; i++) {
    lNorm += evLanczos->data[i] * evLanczos->data[i];
  }

  GenmapGop(c, &lNorm, 1, GENMAP_SCALAR, GENMAP_SUM);
  GenmapScaleVector(evLanczos, evLanczos, 1. / sqrt(lNorm));
  GenmapElements elements = GenmapGetElements(h);
  for(i = 0; i < lelt; i++) {
    elements[i].fiedler = evLanczos->data[i];
  }

  GenmapDestroyVector(alphaVec);
  GenmapDestroyVector(betaVec);
  GenmapDestroyVector(evLanczos);
  GenmapDestroyVector(evTriDiag);
#if defined(GENMAP_TQLI)
  GenmapDestroyVector(eValues);
  for(i = 0; i < iter; i++) {
    GenmapDestroyVector(eVectors[i]);
  }
  GenmapFree(eVectors);
#endif

#if defined(GENMAP_TQLI)
  for(i = 0; i < iter + 1; i++) {
    GenmapDestroyVector(q[i]);
  }
#else
  for(i = 0; i < iter; i++) {
    GenmapDestroyVector(q[i]);
  }
#endif
  GenmapFree(q);

  return iter;
}

int GenmapFiedlerPowerIter(GenmapHandle h,GenmapComm c,
    GenmapVector u,int maxIter, int global)
{
  GenmapVector v;
  GenmapInt lelt=GenmapGetNLocalElements(h);
  GenmapCreateZerosVector(&v,lelt);

  int iter=0;
  GenmapScalar dot;
  dot=GenmapDotVector(v,u);
  GenmapGop(c,&dot,1,GENMAP_SCALAR,GENMAP_SUM);

  GenmapVector weights;
  GenmapCreateVector(&weights,lelt);
  GenmapInitLaplacianWeighted(h,c,weights);

  while(fabs(dot-1)>1e-8 && iter<maxIter){
    GenmapCopyVector(v,u);
    // u=B_g*v
    GenmapLaplacianWeighted(h,c,v,weights,u);

    dot=GenmapDotVector(u,u);
    GenmapGop(c,&dot,1,GENMAP_SCALAR,GENMAP_SUM);
    GenmapScaleVector(u,u,1.0/sqrt(dot));

    dot=GenmapDotVector(v,u);
    GenmapGop(c,&dot,1,GENMAP_SCALAR,GENMAP_SUM);
    iter++;
#if 0
    if(GenmapCommRank(c)==0)
      printf("iter=%03d dot=%lf\n",iter,dot);
#endif
  }

  GenmapInt i;
  GenmapElements e=GenmapGetElements(h);
  for(i=0;i<lelt;i++){
    e[i].fiedler=u->data[i];
  }

  GenmapDestroyVector(weights);
  GenmapDestroyVector(v);

  return iter;
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
      iter = GenmapFiedlerLanczos  (h,c,initVec,maxIter,global);
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
    if(id < (np + 1) / 2) bin = 0;
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
