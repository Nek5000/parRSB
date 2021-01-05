#include <math.h>

#include <genmap-impl.h>
//
//TODO: use a separate function to generate init vector
//
int GenmapFiedlerRQI(genmap_handle h,genmap_comm c,int max_iter,int global)
{
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapVector initVec;
  GenmapCreateVector(&initVec,lelt);

  GenmapElements elements = GenmapGetElements(h);
  GenmapInt i;
  if (global > 0) {
    if (h->options->rsb_paul == 1) {
      for (i = 0;  i < lelt; i++)
        initVec->data[i] = GenmapGetLocalStartIndex(h) + i + 1;
    } else {
      for (i = 0;  i < lelt; i++)
        initVec->data[i] = elements[i].globalId;
    }
  } else {
    for (i = 0;  i < lelt; i++)
      initVec->data[i] = elements[i].fiedler;
  }

  GenmapOrthogonalizebyOneVector(c, initVec, genmap_get_global_nel(h));
  GenmapScalar norm = GenmapDotVector(initVec, initVec);
  GenmapGop(c, &norm, 1, GENMAP_SCALAR, GENMAP_SUM);
  GenmapScalar normi = 1.0/sqrt(norm);
  GenmapScaleVector(initVec, initVec, normi);

  struct comm *gsc=&c->gsc;

  metric_tic(gsc, LAPLACIANSETUP);
  GenmapInitLaplacian(h, c);
  metric_toc(gsc, LAPLACIANSETUP);

  metric_tic(gsc, PRECONDSETUP);
  mgData d;
  mgSetup(h, c, c->M, &d);
  metric_toc(gsc,PRECONDSETUP);

  GenmapVector y;
  GenmapCreateZerosVector(&y,lelt);

  metric_tic(gsc, RQI);
  int iter = rqi(h, c, d, initVec, max_iter, y);
  metric_toc(gsc, RQI);
  metric_acc(NRQI, iter);

  mgFree(d);

  GenmapScalar lNorm = 0;
  for (i = 0; i < lelt; i++)
    lNorm += y->data[i]*y->data[i];
  GenmapGop(c, &lNorm, 1, GENMAP_SCALAR, GENMAP_SUM);

  GenmapScaleVector(y, y, 1./sqrt(lNorm));
  for (i = 0; i < lelt; i++)
    elements[i].fiedler = y->data[i];

  GenmapDestroyVector(y);
  GenmapDestroyVector(initVec);

  return iter;
}

int GenmapFiedlerLanczos(genmap_handle h, genmap_comm c, int max_iter, int global)
{
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapVector initVec, alphaVec, betaVec;

  GenmapCreateVector(&initVec, GenmapGetNLocalElements(h));
  GenmapElements elements = GenmapGetElements(h);

  GenmapInt i;
  if (global > 0) {
    if (h->options->rsb_paul == 1) {
      for(i = 0;  i < lelt; i++)
        initVec->data[i] = GenmapGetLocalStartIndex(h) + i + 1;
    } else {
      for(i = 0;  i < lelt; i++)
        initVec->data[i] = elements[i].globalId;
    }
  } else {
    for(i = 0;  i < lelt; i++)
      initVec->data[i] = elements[i].fiedler;
  }

  GenmapCreateVector(&alphaVec,max_iter);
  GenmapCreateVector(&betaVec,max_iter-1);
  GenmapVector *q = NULL;

  GenmapOrthogonalizebyOneVector(c,initVec,genmap_get_global_nel(h));
  GenmapScalar rtr = GenmapDotVector(initVec, initVec);
  GenmapGop(c, &rtr, 1, GENMAP_SCALAR, GENMAP_SUM);
  GenmapScalar rni = 1.0 / sqrt(rtr);
  GenmapScaleVector(initVec, initVec, rni);

  int iter;
  struct comm *lc = &c->gsc;
  metric_tic(lc, LANCZOS);
  if (h->options->rsb_paul == 1)
    iter = GenmapLanczosLegendary(h, c, initVec, max_iter, &q, alphaVec, betaVec);
  else
    iter = GenmapLanczos(h, c, initVec, max_iter, &q, alphaVec, betaVec);
  metric_toc(lc, LANCZOS);
  metric_acc(NLANCZOS, iter);

  GenmapVector evLanczos, evTriDiag;
  GenmapCreateVector(&evTriDiag, iter);

  /* Use TQLI and find the minimum eigenvalue and associated vector */
  GenmapVector *eVectors, eValues;
  metric_tic(lc, TQLI);
  GenmapTQLI(h, alphaVec, betaVec, &eVectors, &eValues);
  metric_toc(lc, TQLI);

  GenmapScalar eValMin = fabs(eValues->data[0]);
  GenmapInt eValMinI = 0;
  for (i = 1; i < iter; i++) {
    if (fabs(eValues->data[i]) < eValMin) {
      eValMin = fabs(eValues->data[i]);
      eValMinI = i;
    }
  }
  GenmapCopyVector(evTriDiag, eVectors[eValMinI]);

  GenmapInt j;
  GenmapCreateZerosVector(&evLanczos, lelt);
  for (i = 0; i < lelt; i++) {
    for (j = 0; j < iter; j++)
      evLanczos->data[i] += q[j]->data[i] * evTriDiag->data[j];
  }

  GenmapScalar lNorm = 0;
  for (i = 0; i < lelt; i++)
    lNorm += evLanczos->data[i] * evLanczos->data[i];

  GenmapGop(c, &lNorm, 1, GENMAP_SCALAR, GENMAP_SUM);
  GenmapScaleVector(evLanczos, evLanczos, 1. / sqrt(lNorm));
  for(i = 0; i < lelt; i++)
    elements[i].fiedler = evLanczos->data[i];

  GenmapDestroyVector(initVec);
  GenmapDestroyVector(alphaVec);
  GenmapDestroyVector(betaVec);
  GenmapDestroyVector(evLanczos);
  GenmapDestroyVector(evTriDiag);

  if (h->options->rsb_paul == 1) {
    GenmapDestroyVector(eValues);
    for(i = 0; i < iter; i++)
      GenmapDestroyVector(eVectors[i]);
    GenmapFree(eVectors);
  }

  for (i = 0; i < iter; i++)
    GenmapDestroyVector(q[i]);
  if (h->options->rsb_paul == 1)
    GenmapDestroyVector(q[iter]);
  GenmapFree(q);

  return iter;
}

