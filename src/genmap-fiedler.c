#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

#include <genmap-impl.h>

int GenmapFiedler(GenmapHandle h, GenmapComm c, int maxIter, int global){
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapVector initVec, alphaVec, betaVec;

  GenmapCreateVector(&initVec, GenmapGetNLocalElements(h));
  GenmapElements elements = GenmapGetElements(h);

  GenmapInt i;
#if defined(GENMAP_PAUL)
  if(global>0){
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = GenmapGetLocalStartIndex(h) + i + 1;
    }
  }else{
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = elements[i].fiedler;
    }
  }
#else
  if(global>0){
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = (GenmapScalar) elements[i].globalId;
    }
  }else{
    for(i = 0;  i < lelt; i++) {
      initVec->data[i] = elements[i].fiedler;
    }
  }
#endif

  GenmapCreateVector(&alphaVec, maxIter);
  GenmapCreateVector(&betaVec, maxIter - 1);
  GenmapVector *q = NULL;

  GenmapOrthogonalizebyOneVector(h, c, initVec, GenmapGetNGlobalElements(h));
  GenmapScalar rtr = GenmapDotVector(initVec, initVec);
  GenmapGop(c, &rtr, 1, GENMAP_SCALAR, GENMAP_SUM);
  GenmapScalar rni = 1.0 / sqrt(rtr);
  GenmapScaleVector(initVec, initVec, rni);
#if defined(GENMAP_PAUL)
  int iter = GenmapLanczosLegendary(h, c, initVec, maxIter, &q, alphaVec,
                                    betaVec);
#else
  int iter = GenmapLanczos(h, c, initVec, maxIter, &q, alphaVec, betaVec);
#endif
  GenmapVector evLanczos, evTriDiag;
  GenmapCreateVector(&evTriDiag, iter);

#if defined(GENMAP_PAUL)
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

  GenmapInt j;
  GenmapCreateZerosVector(&evLanczos, lelt);
  for(i = 0; i < lelt; i++) {
    for(j = 0; j < iter; j++) {
      evLanczos->data[i] += q[j]->data[i] * evTriDiag->data[j];
    }
  }

  GenmapScalar lNorm = 0;
  for(i = 0; i < lelt; i++) {
    lNorm += evLanczos->data[i] * evLanczos->data[i];
  }

  GenmapGop(c, &lNorm, 1, GENMAP_SCALAR, GENMAP_SUM);
  GenmapScaleVector(evLanczos, evLanczos, 1. / sqrt(lNorm));
  for(i = 0; i < lelt; i++) {
    elements[i].fiedler = evLanczos->data[i];
  }

  GenmapDestroyVector(initVec);
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
