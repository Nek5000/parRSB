#include "genmap-impl.h"

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <stdlib.h>

GenmapScalar g_min,g_max,g_delta;
GenmapLong l_min,l_max;

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

void parRSBGlobalIdMinMax(GenmapHandle h,GenmapComm c,GenmapLong *min,GenmapLong *max) {
  *min = LONG_MAX; *max = LONG_MIN;

  GenmapElements e = GenmapGetElements(h);
  GenmapInt i;
  for(i = 0; i < GenmapGetNLocalElements(h); i++) {
    if(e[i].globalId0 < *min) {
      *min = e[i].globalId0;
    }
    if(e[i].globalId0 > *max) {
      *max = e[i].globalId0;
    }
  }

  GenmapGop(c,min,1,GENMAP_LONG,GENMAP_MIN);
  GenmapGop(c,max,1,GENMAP_LONG,GENMAP_MAX);
}

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

void parRSBHistoSortInitProbes(GenmapHandle h,GenmapComm c,int field) {
  GenmapElements elements=GenmapGetElements(h);
  GenmapInt lelt=GenmapGetNLocalElements(h);

  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);
  int nsplitters=3*(size-1);

  // Allocate space for probes and counts
  GenmapMalloc(nsplitters,&h->histogram->probes);
  GenmapMalloc(nsplitters,&h->histogram->count);

  if(field==GENMAP_FIEDLER)
    parRSBFiedlerMinMax(h,c,&g_min,&g_max);
  else if(field==GENMAP_GLOBALID) {
    parRSBGlobalIdMinMax(h,c,&l_min,&l_max);
    g_max=l_max; g_min=l_min;
  }

  g_delta=(g_max-g_min)/size;
  int i;
  for(i=1; i<size; i++) {
    h->histogram->probes[3*i-3]=g_min;
    h->histogram->probes[3*i-2]=g_min+i*g_delta;
    h->histogram->probes[3*i-1]=g_max;
  }
}

void parRSBHistoSortUpdateCounts(GenmapHandle h,int nsplitters,int field) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  int i;
  for(i=0; i<nsplitters; i++) {
    h->histogram->count[i]=0;
  }

  GenmapElements p,e;
  if(field==GENMAP_FIEDLER){
    for(p=elements,e=p+lelt; p!=e; p++){
      // need to update as we don't keep the probes sorted.
      for(i=0; i<nsplitters; i++){
        if(p->fiedler<h->histogram->probes[i]){
          h->histogram->count[i]++;
        }
      }
    }
  } else if(field==GENMAP_GLOBALID){
    for(p=elements,e=p+lelt; p!=e; p++){
      // need to update as we don't keep the probes sorted.
      for(i=0; i<nsplitters; i++){
        if(p->globalId0<h->histogram->probes[i]){
          h->histogram->count[i]++;
        }
      }
    }
  }
}

int parRSBHistoSortReachedThreshold(GenmapHandle h,GenmapComm c,GenmapLong *count,
                                    GenmapInt threshold)
{
  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);

  GenmapLong lelgt=GenmapGetNGlobalElements(h);
  GenmapInt partition_size=lelgt/size;
  GenmapInt nrem=lelgt-partition_size*size;

  GenmapInt converged=1;

  int i;
  if(rank==0) {
    for(i=1; i<size; i++) {
      GenmapLong expected=i*partition_size+((i<nrem)?i:nrem);
      if(abs(count[3*i-2]-expected)>threshold) {
        converged=0;
        break;
      }
    }
  }

  GenmapBcast(c,&converged,1,GENMAP_INT);
  return converged;
}

void parRSBHistoSortUpdateSplitter(GenmapHandle h,GenmapComm c,int i,int field,
                                   GenmapInt threshold)
{
  int size=GenmapCommSize(c);

  GenmapLong lelgt=GenmapGetNGlobalElements(h);
  GenmapInt partition_size=lelgt/size;
  GenmapInt nrem = lelgt-partition_size*size;

  GenmapLong expected=i*partition_size+((i<nrem)?i:nrem);

  int indx=3*i-2;

  GenmapLong *count=h->histogram->count;
  GenmapScalar *probes=h->histogram->probes;

  if(abs(count[indx]-expected)<threshold) {
    //splitter is achieved
    return;
  }

  if(count[indx]<expected) { // we are going to dump the smaller probe
    probes[indx-1]=probes[indx];
  } else { // we are going to dump the larger probe
    probes[indx+1]=probes[indx];
  }

  probes[indx]=probes[indx-1]+(probes[indx+1]-probes[indx-1])/2;
}

void parRSBHistoSortUpdateProbes(GenmapHandle h,GenmapComm c,
                                 GenmapLong *count,GenmapInt threshold,int field) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  int rank=GenmapCommRank(c);
  int size=GenmapCommSize(c);
  int nsplitters=3*(size-1);

  GenmapLong lelgt = GenmapGetNGlobalElements(h);
  GenmapInt partition_size=lelgt/size;

  int i;
  if(rank==0) {
    for(i=0; i<nsplitters; i++) {
      h->histogram->count[i]=count[i];
    }

    for(i=1; i<size; i++) {
      parRSBHistoSortUpdateSplitter(h,c,i,field,threshold);
    }
  }
}

int parRSBHistoSortSetProc(GenmapHandle h,GenmapComm c,int field,buffer *buf0) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);

  GenmapElements p,e;
  int i;

  if(field==GENMAP_FIEDLER) {
    for(p=elements,e=elements+lelt; p!=e; p++){
      i=1;
      while(i<size && p->fiedler>h->histogram->probes[3*i-2]) i++;
      p->proc=i-1;
    }
  } else if(field==GENMAP_GLOBALID) {
    for(p=elements,e=elements+lelt; p!=e; p++){
      i=1;
      while(i<size && p->globalId0>h->histogram->probes[3*i-2]) i++;
      p->proc=i-1;
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

void parRSBHistoSortLoadBalance(GenmapHandle h,GenmapComm c)
{
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  GenmapLong out[2][1],buf[2][1],in[1];
  GenmapLong n1=lelt;
  in[0]=n1;
  comm_scan(out,&(c->gsComm),genmap_gs_long,gs_add,in,1,buf);

  GenmapLong start=out[0][0];
  GenmapLong nel  =out[1][0];

  GenmapInt id = GenmapCommRank(GenmapGetLocalComm(h));
  GenmapInt np = GenmapCommSize(GenmapGetLocalComm(h));

  GenmapInt pNel = (GenmapInt)(nel / np);
  GenmapInt nrem = (GenmapInt)(nel - pNel * np);
  GenmapInt idCount = 0;
  while(idCount * pNel + ((idCount < nrem) ? idCount : nrem) < start)
    idCount++;

  GenmapLong upLimit = idCount * pNel + ((idCount < nrem) ? idCount : nrem);
  GenmapLong downLimit = start;
  do {
    GenmapInt end = upLimit - start < lelt ? (GenmapInt)(upLimit - start) : lelt;
    GenmapInt i;
    for(i = (GenmapInt)(downLimit - start); i < end; i++)
      elements[i].proc = idCount - 1;
    downLimit = upLimit;
    idCount++;
    upLimit = idCount * pNel + ((idCount < nrem) ? idCount : nrem);
  } while(downLimit - start < lelt);
}

void parRSBHistogramSort(GenmapHandle h,GenmapComm c,int field,buffer *buf0) {
  GenmapScan(h,c);

  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);
  GenmapLong nelg=GenmapGetNGlobalElements(h);

  int print_rank=GenmapCommRank(GenmapGetGlobalComm(h));

  // 10% of load balanced partition size
  GenmapInt threshold=(GenmapGetNGlobalElements(h)/(10*size));
  if(threshold<2) threshold=2;

  // sort locally.
  parRSBHistoSortLocalSort(h,c,field,buf0);

  // We are done if size==1
  if(size==1) return;

  // Else we continue
  int nsplitters=3*(size-1);

  // Allocate space for probes and counts
  GenmapMalloc(nsplitters,&(h->histogram->probes));
  GenmapMalloc(nsplitters,&(h->histogram->count));
  GenmapLong *count=NULL;
  if(rank==0) {
    GenmapMalloc(nsplitters,&count);
  }

  // init probes values
  parRSBHistoSortInitProbes(h,c,field);
#if 0
  if(print_rank==0)
    for(int i=0; i<nsplitters; i++) {
      printf("iter: 0 probe[%d]= " GenmapScalarFormat "\n",i,h->histogram->probes[i]);
    }
#endif

  // update counts locally
  parRSBHistoSortUpdateCounts(h,nsplitters,field);

  // global reduction
  GenmapReduce(c,count,h->histogram->count,nsplitters,GENMAP_LONG,GENMAP_SUM);
#if 0
  if(print_rank==0)
    for(int i=0; i<nsplitters; i++) {
      printf("iter: 0 count[%d]= " GenmapLongFormat "\n",i,count[i]);
    }
#endif

  int iter=0;
  while(!parRSBHistoSortReachedThreshold(h,c,count,threshold)){
    parRSBHistoSortUpdateProbes(h,c,count,threshold,field);
#if 0
    if(print_rank==0)
      for(int i=1; i<size; i++){
        printf("iter: %d probe[%d]=" GenmapScalarFormat " count[%d]=" GenmapLongFormat
            "\n",iter+1,i,h->histogram->probes[3*i-2],i,h->histogram->count[3*i-2]);
      }
#endif

    // TODO: Bcast probes
    GenmapBcast(c,h->histogram->probes,nsplitters,GENMAP_LONG);
#if 0
    if(print_rank==0)
      for(int i=0; i<nsplitters; i++){
        printf("%d: %d probe[%d]= " GenmapScalarFormat "\n",print_rank,iter,i,
            h->histogram->probes[i]);
      }
#endif

    parRSBHistoSortUpdateCounts(h,nsplitters,field);
#if 0
    if(print_rank==0)
      for(int i=1; i<size; i++){
        printf("iter: %d probe[%d]=" GenmapScalarFormat " count[%d]=" GenmapLongFormat
            "\n",iter+1,i,h->histogram->probes[3*i-2],i,h->histogram->count[3*i-2]);
      }
#endif

    // global reduction
    GenmapReduce(c,count,h->histogram->count,nsplitters,GENMAP_LONG,GENMAP_SUM);
#if 0
    if(print_rank==0)
      for(int i=1; i<size; i++){
        printf("iter: %d probe[%d]=" GenmapScalarFormat " count[%d]=" GenmapLongFormat
            " min=%lf max=%lf size=%d nelg=%lld\n",iter+1,i,h->histogram->probes[3*i-2],i,
            h->histogram->count[3*i-2],g_min,g_max,size,nelg);
      }
#endif
    iter++;
#if 0
    if(iter>=10) {
      MPI_Finalize();
      exit(0);
    }
#endif
  }
  // set destination processor id for each element
  parRSBHistoSortSetProc(h,c,field,buf0);

  // send elements to right processor
  parRSBHistoSortTransferToProc(h,field,buf0);
  GenmapScan(h,c);
  // sort locally.
  parRSBHistoSortLocalSort(h,c,field,buf0);

  GenmapInt lelt = GenmapGetNLocalElements(h);
  GenmapInt zero=0,sum=0;
  if(lelt<=0) {
    zero=1;
    sum=1;
  }
  GenmapGop(c,&sum,1,GENMAP_INT,GENMAP_SUM);
  GenmapGop(c,&lelt,1,GENMAP_INT,GENMAP_SUM);
  if(print_rank==0) printf("sum0=%d nelg=%d size=%d\n",sum,lelt,size);

  // do a load balancing step
  parRSBHistoSortLoadBalance(h,c);

  // send elements to right processor
  parRSBHistoSortTransferToProc(h,field,buf0);
  GenmapScan(h,c);
  parRSBHistoSortLocalSort(h,c,field,buf0);

  lelt = GenmapGetNLocalElements(h);
  zero=0,sum=0;
  if(lelt<=0) {
    zero=1;
    sum=1;
  }
  GenmapGop(c,&sum,1,GENMAP_INT,GENMAP_SUM);
  GenmapGop(c,&lelt,1,GENMAP_INT,GENMAP_SUM);
  if(print_rank==0) printf("sum1=%d nelg=%d size=%d\n",sum,lelt,size);

  // Finalize sort
  GenmapFree(h->histogram->probes);
  GenmapFree(h->histogram->count);
  if(rank==0) {
    GenmapFree(count);
  }
#if 0
    if(print_rank==0)
      printf("Done with histo-sort\n");
#endif
}