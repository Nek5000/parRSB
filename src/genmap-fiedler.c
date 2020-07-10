#include "genmap-impl.h"

#include <math.h>
#include <stdio.h>

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
