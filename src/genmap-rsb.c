#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

#include <genmap-impl.h>
#include <sort.h>
#include <parRSB.h>

void genmap_rsb(genmap_handle h,int verbose){
  int max_iter=50;
  int max_pass=50;

  GenmapComm local_c=GenmapGetLocalComm(h);
  struct comm *lc=&local_c->gsc;

  GenmapComm global_c=GenmapGetGlobalComm(h);
  struct comm *gc=&global_c->gsc;

  GenmapScan(h,local_c);
  crystal_init(&h->cr,lc);

  GenmapScan(h, GenmapGetLocalComm(h));
  uint nelt=GenmapGetNLocalElements(h);
  GenmapElements e=GenmapGetElements(h);
  GenmapInt i;
  for(i=0; i<nelt; i++) {
    e[i].globalId0=GenmapGetLocalStartIndex(h)+i+1;
  }

  buffer buf0=null_buffer;

  int nve=h->nv;
  int ndim=(nve==8)?3:2;
  int level=0;
  int np;

  while(GenmapCommSize(GenmapGetLocalComm(h))>1){
    local_c=GenmapGetLocalComm(h);
    lc=&local_c->gsc;
    np=lc->np;

    metric_tic(lc,RSB);

#if defined(GENMAP_PAUL)
    int global=1;
#else
    int global=(np==gc->np);
#endif

    /* Run RCB pre-step */
    rcb(lc,h->elementArray,ndim);

    /* Initialize the laplacian */
    nelt=GenmapGetNLocalElements(h);
    GenmapCreateVector(&h->weights,nelt);
    GenmapInitLaplacianWeighted(h,local_c,h->weights);

    /* Run fiedler */
    metric_tic(lc,FIEDLER);

    int ipass=0,iter;
    do {
#if defined(GENMAP_LANCZOS)
      iter=GenmapFiedlerLanczos(h,local_c,max_iter,global);
#elif defined(GENMAP_RQI)
      iter=GenmapFiedlerRQI(h,local_c,max_iter,global);
#endif
      metric_acc(NFIEDLER,iter);
      global=0;
    } while(++ipass<max_pass && iter==max_iter);
    metric_toc(lc,FIEDLER);

    parallel_sort(struct rsb_element,h->elementArray,fiedler,
      gs_double,0,1,lc);

    metric_toc(lc,RSB);

    GenmapInt id=lc->id;
    int bin=1; if(id<(np+1)/2) bin=0;
    GenmapSplitComm(h,&local_c,bin);
    GenmapSetLocalComm(h,local_c);
    lc=&local_c->gsc;

    GenmapScan(h,local_c);

    // FIXME: Do this only at begining of the loop
    GenmapFree(h->weights);

    level++;
    metric_push_level();
  }

  /* Check if Fidler converged */
  sint converged=1,buf;
  for(i=0; i<metric_get_levels(); i++){
    int val=(int)metric_get_value(i,NFIEDLER);
    if(val>=max_pass*max_iter){
      converged=0; break;
    }
  }
  comm_allreduce(gc,gs_int,gs_min,&converged,1,&buf);// min
  if(converged==0 && gc->id==0)
    printf("\tWARNING: Lanczos failed to converge while partitioning!\n");

  /* Check for disconnected components */
  nelt=GenmapGetNLocalElements(h);
  GenmapCreateVector(&h->weights,nelt);
  GenmapInitLaplacianWeighted(h,global_c,h->weights);

  sint discon=is_disconnected(gc,global_c->gsh,&global_c->buf,nelt,nve);
  if(discon>0 && gc->id==0)
    printf("\tWarning: There are disconnected components!\n");

  GenmapFree(h->weights);

  crystal_free(&h->cr);
  buffer_free(&buf0);
}
