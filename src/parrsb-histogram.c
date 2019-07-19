#include "genmap-impl.h"

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

void parRSBHistoSortInitProbes(GenmapHandle h,int field) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  GenmapComm c=GenmapGetGlobalComm(h);

  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);
  int nsplitters=3*(size-1);

  if(field==GENMAP_FIEDLER){
    // pick 3*(size-1) splitters as probes
    GenmapScalar delta=(1.f-(-1.f))/size;
    for(int i=1; i<size; i++) {
      h->histogram->probes[3*i-3]=-1.f+i*delta-delta/4;
      h->histogram->probes[3*i-2]=-1.f+i*delta;
      h->histogram->probes[3*i-1]=-1.f+i*delta+delta/4;
    }

    // initialize the count to zero
    for(int i=0; i<nsplitters; i++) {
      h->histogram->count[i] = 0;
    }
  }
}

void parRSBHistoSortUpdateCounts(GenmapHandle h,int field) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  GenmapComm c=GenmapGetGlobalComm(h);
  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);
  int nsplitters=3*(size-1);

  if(field == GENMAP_FIEDLER) {
    // calculate initial counts
    // local
    GenmapElements p,e;
    for(p=elements,e=elements+lelt; p!=e; p++) {
      for(int i=0; i<nsplitters; i++) {
        if(p->fiedler<h->histogram->probes[i]) {
          h->histogram->count[i]++;
        }
      }
    }
  }

}

void parRSBHistoSortLocalSort(GenmapHandle h,int field,buffer *buf0) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  GenmapComm c=GenmapGetGlobalComm(h);
  if(field == GENMAP_FIEDLER) { // Fiedler
    // sort locally according to Fiedler vector
    sarray_sort(struct GenmapElement_private, elements, (GenmapUInt)lelt, fiedler,
                TYPE_DOUBLE, buf0);
  } else if(GENMAP_GLOBALID) {
    // sort locally according to globalId
    sarray_sort(struct GenmapElement_private, elements, (GenmapUInt)lelt,
                globalId0, TYPE_LONG, buf0);
  }

}

int parRSBHistoSortReachedThreshold(GenmapHandle h,GenmapInt threshold,int field) {
  GenmapComm c=GenmapGetGlobalComm(h);
  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);

  GenmapLong lelgt = GenmapGetNGlobalElements(h);
  GenmapInt partition_size=lelgt/size;

  GenmapLong *current = h->histogram->count;
  GenmapScalar *probes = h->histogram->probes;

  GenmapInt converged=1;

  if(rank==0) {
    for(int i=1; i<size; i++) {
      if(abs(current[3*i-2]-partition_size)>threshold) {
        converged=0;
        break;
      }
    }
  }

  GenmapBcast(c,&converged,1,GENMAP_INT);
  return converged;
}

void parRSBHistoSortUpdateProbes(GenmapHandle h,GenmapLong *count,GenmapInt threshold,
                                 int field) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  GenmapComm c=GenmapGetGlobalComm(h);
  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);
  int nsplitters=3*(size-1);

  GenmapLong lelgt = GenmapGetNGlobalElements(h);
  GenmapInt partition_size=lelgt/size;

  GenmapLong *current = h->histogram->count;
  GenmapScalar *probes = h->histogram->probes;

  GenmapScalar gradient;

  if(rank==0) {
    if(field == GENMAP_FIEDLER) {
      for(int i=1; i<size; i++) {
        if(abs(current[3*i-2]-partition_size)>abs(count[3*i-2]-partition_size)) {
          //update
          current[3*i-3]=count[3*i-3];
          current[3*i-2]=count[3*i-2];
          current[3*i-1]=count[3*i-1];
          // adjust middle probe
          gradient=(probes[3*i-1]-probes[3*i-3])/(count[3*i-1]-count[3*i-3]);
          probes[3*i-2]=probes[3*i-3]+(partition_size-count[3*i-3])*gradient;
        }
        // adjust left probe
        if((current[3*i-2]-current[3*i-3])>threshold/2) {
          gradient=(probes[3*i-2]-probes[3*i-3])/(current[3*i-2]-current[3*i-3]);
          probes[3*i-3]=probes[3*i-2]-(threshold/2)*gradient;
        }
        // adjust right probe
        if((current[3*i-1]-current[3*i-2])>threshold/2) {
          gradient=(probes[3*i-1]-probes[3*i-2])/(current[3*i-1]-current[3*i-2]);
          probes[3*i-1]=probes[3*i-2]+(threshold/2)*gradient;
        }
      }
    }
  }
}


void parRSBHistogramSort(GenmapHandle h,GenmapComm c,int field,buffer *buf0) {
  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);
  int nsplitters=3*(size-1);

  // Allocate space for probes and counts
  GenmapMalloc(nsplitters,&h->histogram->probes);
  GenmapMalloc(nsplitters,&h->histogram->count);
  GenmapLong *count;
  if(rank==0) {
    GenmapMalloc(nsplitters,&count);
  }

  // 10% of load balanced partition size 
  GenmapInt threshold=(GenmapGetNGlobalElements(h)/size)/10;
  if(threshold<2) threshold=2;

  // sort locally.
  parRSBHistoSortLocalSort(h,field,buf0);
  // init probes values
  parRSBHistoSortInitProbes(h,field);
  // update counts locally 
  parRSBHistoSortUpdateCounts(h,field);

  // global reduction
  GenmapReduce(c,count,h->histogram->count,nsplitters,GENMAP_LONG,GENMAP_SUM);

  while(!parRSBHistoSortReachedThreshold(h,threshold,field)) {
     parRSBHistoSortUpdateProbes(h,count,threshold,field);

     // TODO: Bcast probes
     GenmapBcast(c,h->histogram->probes,nsplitters,GENMAP_LONG);

     parRSBHistoSortUpdateCounts(h,field);

     // global reduction
     GenmapReduce(c,count,h->histogram->count,nsplitters,GENMAP_LONG,GENMAP_SUM);
  }

  GenmapFree(h->histogram->probes);
  GenmapFree(h->histogram->count);
  if(rank==0) {
    GenmapFree(count);
  }
}
