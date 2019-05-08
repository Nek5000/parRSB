#include <genmap-impl.h>

int parRSBInitLaplacian(GenmapHandle h) {
  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapInt nv = GenmapGetNVertices(h);
  GenmapUInt numPoints = (GenmapUInt) nv * lelt;

  GenmapLong *vertices;
  GenmapMalloc(numPoints, &vertices);

  GenmapElements elements = GenmapGetElements(h);
  GenmapInt i, j;
  for(i = 0; i < lelt; i++) {
    for(j = 0; j < nv; j++) {
      vertices[i * nv + j] = elements[i].vertices[j];
    }
  }

  if(h->local->verticesHandle)
    gs_free(h->local->verticesHandle);

#if defined(GENMAP_DEBUG)
  double t1 = GenmapGetMaxRss();
  if(GenmapCommRank(GenmapGetLocalComm(h)) == 0)
    printf("RSS before gs_setup: %lf\n", t1);
#endif

  h->local->verticesHandle = gs_setup(vertices, numPoints, &h->local->gsComm, 0,
                                      gs_crystal_router, 0);
#if defined(GENMAP_DEBUG)
  t1 = GenmapGetMaxRss();
  if(GenmapCommRank(GenmapGetLocalComm(h)) == 0)
    printf("RSS after gs_setup: %lf\n", t1);
#endif

  GenmapScalar *u;
  GenmapMalloc(numPoints, &u);

  for(i = 0; i < lelt; i++)
    for(j = 0; j < nv; j++)
      u[nv * i + j] = 1.;

  gs(u, genmap_gs_scalar, gs_add, 0, h->local->verticesHandle, &h->local->buf);

  GenmapCreateVector(&h->weights, lelt);

  for(i = 0; i < lelt; i++) {
    h->weights->data[i] = 0.;
    for(j = 0; j < nv; j++) {
      h->weights->data[i] += u[nv * i + j];
    }
  }

  for(i = 0; i < lelt; i++) {
    h->weights->data[i] *= -1;
  }

  GenmapFree(u);
  GenmapFree(vertices);

  return 0;
}

int parRSBLaplacian(GenmapHandle h, GenmapVector u, GenmapVector v) {
  assert(u->size == v->size);
  assert(u->size == GenmapGetNLocalElements(h));

  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapInt nv = GenmapGetNVertices(h);

  GenmapScalar *ucv;
  GenmapMalloc((size_t)(nv * lelt), &ucv);

  GenmapInt i, j;
  for(i = 0; i < lelt; i++)
    for(j = 0; j < nv; j++)
      ucv[nv * i + j] = u->data[i];

  gs(ucv, genmap_gs_scalar, gs_add, 0, h->local->verticesHandle, &h->local->buf);

  for(i = 0; i < lelt; i++) {
    v->data[i] = h->weights->data[i] * u->data[i];
    for(j = 0; j < nv; j ++) {
      v->data[i] += ucv[nv * i + j];
    }
  }

  GenmapFree(ucv);

  return 0;
}

int parRSBFinalizeLaplacian(GenmapHandle h) {
  GenmapDestroyVector(h->weights);
}
