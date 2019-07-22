#include "genmap-impl.h"

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

// A comparator function used by qsort 
int compare(const void * a, const void * b) 
{ 
    return ( *(GenmapScalar*)a - *(GenmapScalar*)b ); 
} 

GenmapScalar g_min,g_max;

void parRSBHistoSortLocalSort(GenmapHandle h,GenmapComm c,int field,buffer *buf0) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

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

int parRSBHistoSortReachedThreshold(GenmapHandle h,GenmapComm c,GenmapLong *count,
                                    GenmapInt threshold,int field) {
  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);

  GenmapLong lelgt = GenmapGetNGlobalElements(h);
  GenmapInt partition_size=lelgt/size;
  int nsplitters=size-1;

  GenmapInt converged=1;

  if(rank==0) {
    for(int i=0; i<nsplitters; i++) {
      if(abs(count[i]-i*partition_size)>threshold) {
        converged=0;
        break;
      }
    }
  }

  GenmapBcast(c,&converged,1,GENMAP_INT);
  return converged;
}

void parRSBFiedlerMinMax(GenmapHandle h,GenmapComm c,GenmapScalar *min,GenmapScalar *max) {
  *min = 1; *max = -1;

  GenmapElements e = GenmapGetElements(h);
  GenmapInt i;
  for(i = 0; i < GenmapGetNLocalElements(h); i++) {
    if(e[i].fiedler < *min) {
      *min = e[i].fiedler;
    }
    if(e[i].fiedler > *max) {
      *max = e[i].fiedler;
    }
  }

  GenmapGop(c,min,1,GENMAP_SCALAR,GENMAP_MIN);
  GenmapGop(c,max,1,GENMAP_SCALAR,GENMAP_MAX);
}

void parRSBHistoSortInitProbes(GenmapHandle h,GenmapComm c,int field) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);
  int nsplitters=size-1;

  // Allocate space for probes and counts
  GenmapMalloc(nsplitters,&h->histogram->probes);
  GenmapMalloc(nsplitters,&h->histogram->count);

  parRSBFiedlerMinMax(h,c,&g_min,&g_max);

  if(field==GENMAP_FIEDLER){
    // pick 3*(size-1) splitters as probes
    GenmapScalar delta=(g_max-g_min)/size;
    for(int i=0; i<nsplitters; i++) {
      h->histogram->probes[i]=g_min+(i+1)*delta;
    }

    // initialize the count to zero
    for(int i=0; i<nsplitters; i++) {
      h->histogram->count[i] = 0;
    }
  }
}

void parRSBHistoSortUpdateCounts(GenmapHandle h,int nsplitters,int field) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  for(int i=0; i<nsplitters; i++) {
    h->histogram->count[i]=0;
  }

  if(field == GENMAP_FIEDLER) {
    // calculate initial counts
    // local
    GenmapElements p,e;
    for(p=elements,e=p+lelt; p!=e; p++) {
      int i=0;
      while(p->fiedler>h->histogram->probes[i] && i<nsplitters) i++;
      while(p->fiedler<=h->histogram->probes[i] && i<nsplitters) {
        h->histogram->count[i]++;
        i++;
      }
    }
  }
}

void parRSBHistoSortUpdateProbes(GenmapHandle h,GenmapComm c,
                                 GenmapLong *count,GenmapInt threshold,int field) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);
  int nsplitters=size-1;

  GenmapLong lelgt = GenmapGetNGlobalElements(h);
  GenmapInt partition_size=lelgt/size;

  GenmapLong *current = h->histogram->count;
  GenmapScalar *probes = h->histogram->probes;

  GenmapScalar gradient;

  if(rank==0) {
    if(field == GENMAP_FIEDLER) {
      for(int i=0; i<nsplitters; i++) {
        if(abs(current[i]-(i+1)*partition_size)>abs(count[i]-(i+1)*partition_size)) {
          //update
          //printf("I am here: %d %lld %lld",i,current[3*i-2],count[3*i-2]);
          current[i]=count[i];
        }
      }
      // Now adjust probes based on linear spline
      GenmapScalar p_probes=g_min;
      GenmapLong p_count=0;
      for(int i=0; i<nsplitters; i++) {
        GenmapScalar old_probe=h->histogram->probes[i];
        GenmapScalar gradient=(count[i]-p_count)/(probes[i]-p_probes);
        h->histogram->probes[i]=((i+1)*partition_size-p_count)/gradient+p_probes;
        p_count=count[i];
        p_probes=old_probe;
      }
      // Not sort the probes
      qsort(h->histogram->probes,nsplitters,sizeof(GenmapScalar),compare); 
    }
  }
}

int parRSBHistoSortSetProc(GenmapHandle h,GenmapComm c,int field,buffer *buf0) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);
  int nsplitters=size-1;
  
  GenmapElements p,e;
  int i;

  if(field==GENMAP_FIEDLER) {
    for(p=elements,e=elements+lelt; p!=e; p++){
      i=0;
      while(i<nsplitters && p->fiedler>h->histogram->probes[i]) i++;
      p->proc=i;
    }
  }

  return 0;
}

void parRSBHistoSortTransferToProc(GenmapHandle h, int field, buffer *buf0) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  if(field == GENMAP_FIEDLER) { // Fiedler
    sarray_transfer(struct GenmapElement_private, &(h->elementArray), proc, 0,
                    &(h->cr));
    GenmapScan(h, GenmapGetLocalComm(h));
  } else if(field == GENMAP_GLOBALID) {
    sarray_transfer(struct GenmapElement_private, &(h->elementArray), proc, 0,
                    &(h->cr));
    GenmapScan(h, GenmapGetLocalComm(h));
  }
}

void parRSBHistogramSort(GenmapHandle h,GenmapComm c,int field,buffer *buf0) {
  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);

  // 10% of load balanced partition size 
  GenmapInt threshold=(GenmapGetNGlobalElements(h)/size);
  if(threshold<2) threshold=2;

  // sort locally.
  parRSBHistoSortLocalSort(h,c,field,buf0);

  // We are done if size==1
  if(size==1) return;

  // Else we continue
  int nsplitters=size-1;

  // Allocate space for probes and counts
  GenmapMalloc(nsplitters,&(h->histogram->probes));
  GenmapMalloc(nsplitters,&(h->histogram->count));
  GenmapLong *count=NULL;
  if(rank==0) {
    GenmapMalloc(nsplitters,&count);
  }

  // init probes values
  parRSBHistoSortInitProbes(h,c,field);
  if(rank==0)
    for(int i=0; i<nsplitters; i++) {
      printf("%d -1 probe[%d]= " GenmapScalarFormat "\n",rank,i,h->histogram->probes[i]);
    }

  // update counts locally 
  parRSBHistoSortUpdateCounts(h,nsplitters,field);

  // global reduction
  GenmapReduce(c,count,h->histogram->count,nsplitters,GENMAP_LONG,GENMAP_SUM);
  if(rank==0)
    for(int i=0; i<nsplitters; i++) {
      printf("iter: -1 count[%d]= " GenmapLongFormat "\n",i,count[i]);
    }

  int iter=0;
  while(!parRSBHistoSortReachedThreshold(h,c,count,threshold,field)){
    parRSBHistoSortUpdateProbes(h,c,count,threshold,field);
    if(rank==0)
      for(int i=0; i<nsplitters; i++){
        printf("%d: %d probe[%d]= " GenmapScalarFormat "\n",rank,iter,i,h->histogram->probes[i]);
      }
    if(rank==0) printf("Update probes: done\n");

    // TODO: Bcast probes
    GenmapBcast(c,h->histogram->probes,nsplitters,GENMAP_LONG);
    if(rank==1)
      for(int i=0; i<nsplitters; i++){
        printf("%d: %d probe[%d]= " GenmapScalarFormat "\n",rank,iter,i,h->histogram->probes[i]);
      }
    if(rank==0) printf("Bcast: done\n");

    parRSBHistoSortUpdateCounts(h,nsplitters,field);
    if(rank==0)
      for(int i=0; i<nsplitters; i++){
        printf("%d: %d count[%d]= " GenmapLongFormat "\n",rank,iter,i,h->histogram->count[i]);
      }
    if(rank==0) printf("Update counts: done\n");

    // global reduction
    GenmapReduce(c,count,h->histogram->count,nsplitters,GENMAP_LONG,GENMAP_SUM);
    if(rank==0) printf("Global reduce: done\n");
    if(rank==0)
      for(int i=0; i<nsplitters; i++){
        printf("iter: %d count[%d]= " GenmapLongFormat "\n",iter,i,count[i]);
      }
    iter++;
    if(iter>10) {
      MPI_Finalize();
      exit(0);
    }
  }
  // set destination processor id for each element
  parRSBHistoSortSetProc(h,c,field,buf0);

  // send elements to right processor
  GenmapCrystalInit(h,c);
  parRSBHistoSortTransferToProc(h,field,buf0);
  GenmapCrystalFinalize(h);

  GenmapScan(h,c);

  // sort locally.
  parRSBHistoSortLocalSort(h,c,field,buf0);

  // Finalize sort
  GenmapFree(h->histogram->probes);
  GenmapFree(h->histogram->count);
  if(rank==0) {
    GenmapFree(count);
  }
}
