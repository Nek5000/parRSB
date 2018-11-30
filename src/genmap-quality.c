#include <genmap-impl.h>

#include <math.h>
#include <stdio.h>

GenmapInt GenmapPartitionQuality(GenmapHandle h) {
  GenmapInt id = h->Id(h->global);
  GenmapInt np = h->Np(h->global);
  GenmapInt lelt = h->header->lelt;
  GenmapInt nv = h->header->nv;

  GenmapLong *data;
  GenmapInt numPoints = lelt * nv;
  GenmapMalloc(numPoints, &data);

  GenmapElements elements = GenmapGetElements(h);
  for(GenmapInt i = 0; i < lelt; i++) {
    for(int j = 0; j < nv; j++) {
      data[i * nv + j] = elements[i].vertices[j];
    }
  }

  GenmapComm c;
  GenmapCreateComm(&c, h->global->gsComm.c);
  c->verticesHandle = gs_setup(data, numPoints, &c->gsComm, 0, gs_pairwise,
                               0);

  GenmapInt neighborsCount = 0;
  for(GenmapInt i = 0; i < np; i++) {
    if(i != id) {
      for(GenmapInt j = 0; j < numPoints; j++) {
        data[j] = -1;
      }
    } else {
      for(GenmapInt j = 0; j < numPoints; j++) {
        data[j] = id + 1;
      }
    }

    gs(data, gs_long, gs_max, 0, c->verticesHandle, NULL);

    for(GenmapInt j = 0; j < numPoints; j++) {
      if(data[j] > 0) {
        neighborsCount++;
        break;
      }
    }
  }

  GenmapInt ncMax = neighborsCount;
  GenmapInt ncMin = neighborsCount;
  GenmapInt ncSum = neighborsCount;

  GenmapGop(c, &ncMax, 1, GENMAP_INT, GENMAP_MAX);
  GenmapGop(c, &ncMin, 1, GENMAP_INT, GENMAP_MIN);
  GenmapGop(c, &ncSum, 1, GENMAP_INT, GENMAP_SUM);

  if(GenmapId(h->global) == 0) {
    printf("Max neighbors: "GenmapIntFormat,ncMax);
    printf(" Min neighbors: "GenmapIntFormat,ncMin);
    printf(" Avg neighbors: "GenmapScalarFormat"\n",(1.0*ncSum)/GenmapNp(h->global));
  }

  GenmapFree(data);
  GenmapDestroyComm(c); 

  return neighborsCount;
}