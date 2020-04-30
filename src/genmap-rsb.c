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
#if defined(GENMAP_PAUL)
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

#if defined(GENMAP_PAUL)
  int iter = GenmapLanczosLegendary(h, c, initVec, maxIter, &q, alphaVec,
                                    betaVec);
#else
  int iter = GenmapLanczos         (h, c, initVec, maxIter, &q, alphaVec,
                                    betaVec);
#endif
  GenmapVector evLanczos, evTriDiag;
  GenmapCreateVector(&evTriDiag, iter);

#if defined(GENMAP_PAUL)
  /* Use TQLI and find the minimum eigenvalue and associated vector */
  GenmapVector *eVectors, eValues;
  GenmapTQLI(h, alphaVec, betaVec, &eVectors, &eValues);

  GenmapScalar eValMin = fabs(eValues->data[0]);
  GenmapInt eValMinI = 0;
  GenmapInt i;
  for(i = 1; i < iter; i++) {
    if(fabs(eValues->data[i]) < eValMin) {
      eValMin = fabs(eValues->data[i]);
      eValMinI = i;
    }
  }
  GenmapCopyVector(evTriDiag, eVectors[eValMinI]);
#if defined(GENMAP_DEBUG)
  if(GenmapCommRank(GenmapGetGlobalComm(h)) == 0) {
    printf("evTriDiag:\n");
    GenmapPrintVector(evTriDiag);
  }
#endif
#else
  GenmapVector init;
  GenmapCreateVector(&init, iter);
  for(int i = 0; i < iter; i++) {
    init->data[i] = i + 1.0;
  }
  GenmapScalar avg = 0.5 * iter * (1.0 + iter) / iter;
  for(int i = 0; i < iter; i++) {
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
    for(int i = 0; i < evLanczos->size; i++) {
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
#if defined(GENMAP_PAUL)
  GenmapDestroyVector(eValues);
  for(i = 0; i < iter; i++) {
    GenmapDestroyVector(eVectors[i]);
  }
  GenmapFree(eVectors);
#endif

#if defined(GENMAP_PAUL)
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
    GenmapVector initVec,int maxIter, int global)
{
  //return iter;
}

int GenmapFiedler(GenmapHandle h,GenmapComm c,int maxIter,int global){
  GenmapVector initVec;
  GenmapGetInitVector(h,c,global,&initVec);
  int iter=GenmapFiedlerLanczos(h,c,initVec,maxIter,global);
  GenmapDestroyVector(initVec);
  return iter;
}

void GenmapRSB(GenmapHandle h) {
  int maxIter = 50;
  int npass = 50;

  GenmapInt i;
  GenmapElements e = GenmapGetElements(h);
  GenmapScan(h, GenmapGetLocalComm(h));
  for(i = 0; i < GenmapGetNLocalElements(h); i++) {
    e[i].globalId  = GenmapGetLocalStartIndex(h) + i + 1;
    e[i].globalId0 = GenmapGetLocalStartIndex(h) + i + 1;
  }

  if(GenmapCommRank(GenmapGetGlobalComm(h)) == 0 && h->dbgLevel > 0)
    printf("running RSB ...\n"), fflush(stdout);

  crystal_init(&(h->cr), &(h->local->gsComm));
  buffer buf0 = null_buffer;

  int level=0;

  while(GenmapCommSize(GenmapGetLocalComm(h)) > 1) {
#if defined(GENMAP_PAUL)
    int global = 1;
#else
    int global=(GenmapCommSize(GenmapGetLocalComm(h))== \
        GenmapCommSize(GenmapGetGlobalComm(h)));
#endif

    int ipass = 0;
    int iter;
    do {
      iter = GenmapFiedler(h, GenmapGetLocalComm(h), maxIter, global);
      ipass++;
      global = 0;
    } while(ipass < npass && iter == maxIter);

    if(GenmapCommRank(GenmapGetGlobalComm(h))==0 && h->dbgLevel>1){
      printf("level %02d: %5d iterations\n",level,(ipass-1)*maxIter+iter);
      fflush(stdout);
    }

    GenmapBinSort(h, GENMAP_FIEDLER, &buf0);

    int bin;
    GenmapInt np = GenmapCommSize(GenmapGetLocalComm(h));
    GenmapInt id = GenmapCommRank(GenmapGetLocalComm(h));
    if(id < (np + 1) / 2) bin = 0;
    else bin = 1;

    GenmapComm c = GenmapGetLocalComm(h);
    GenmapSplitComm(h, &c, bin);
    GenmapSetLocalComm(h, c);

#if defined(GENMAP_PAUL)
    GenmapBinSort(h, GENMAP_GLOBALID, &buf0);
#endif
    level++;
  }

  crystal_free(&(h->cr));
  buffer_free(&buf0);
}
