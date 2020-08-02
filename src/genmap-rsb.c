#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

#include <genmap-impl.h>

void GenmapRSB(GenmapHandle h,int verbose){
  int maxIter=50;
  int npass  =50;

  GenmapInt i;
  GenmapElements e = GenmapGetElements(h);
  GenmapScan(h, GenmapGetLocalComm(h));
  for(i = 0; i < GenmapGetNLocalElements(h); i++) {
    e[i].globalId =GenmapGetLocalStartIndex(h)+i+1;
    e[i].globalId0=GenmapGetLocalStartIndex(h)+i+1;
  }

  if(GenmapCommRank(GenmapGetGlobalComm(h)) == 0 && h->dbgLevel > 0)
    printf("running RSB "), fflush(stdout);

  crystal_init(&(h->cr), &(h->local->gsc));
  buffer buf0 = null_buffer;

  int rank=GenmapCommRank(GenmapGetGlobalComm(h));
  int level=0;

  while(GenmapCommSize(GenmapGetLocalComm(h)) > 1){
    GenmapComm local_c=GenmapGetLocalComm(h);
    GenmapInt np=GenmapCommSize(local_c);

#if defined(GENMAP_PAUL)
    int global=1;
#else
    int global=(np==GenmapCommSize(GenmapGetGlobalComm(h)));
#endif

    for(i=0; i<16; i++)
      h->time[i]=0.0;

    int ipass=0,iter;
    do{
      comm_barrier(&local_c->gsc);
      h->time[0]-=comm_time();
#if 1
      iter=GenmapFiedlerLanczos(h,local_c,maxIter,global);
#else
      iter=GenmapFiedlerRQI(h,local_c,maxIter,global);
#endif
      comm_barrier(&local_c->gsc);
      h->time[0]+=comm_time();
      h->time[1]+=iter;
      global=0;
    }while(++ipass < npass && iter==maxIter);

    double min[16],max[16],sum[16],buf[16];
    for(i=0; i<16; i++)
      min[i]=max[i]=sum[i]=h->time[i];

    comm_allreduce(&local_c->gsc,gs_double,gs_min,min,16,buf); // min
    comm_allreduce(&local_c->gsc,gs_double,gs_max,max,16,buf); // max
    comm_allreduce(&local_c->gsc,gs_double,gs_add,sum,16,buf); // sum

    if(rank==0 && verbose){
      printf("level=%02d:\n",level);
      printf("fiedler(time)  : %lf/%lf/%lf\n",min[0],max[0],sum[0]/np);
      printf("fiedler(iter)  : %lf/%lf/%lf\n",min[1],max[1],sum[1]/np);
      printf("fine-op setup  : %lf/%lf/%lf\n",min[2],max[2],sum[2]/np);
      printf("fine-op apply  : %lf/%lf/%lf\n",min[11],max[11],sum[11]/np);
      printf("rqi(time)      : %lf/%lf/%lf\n",min[3],max[3],sum[3]/np);
      printf("rqi(iter)      : %lf/%lf/%lf\n",min[12],max[12],sum[12]/np);
      printf("mg setup       : %lf/%lf/%lf\n",min[4],max[4],sum[4]/np);
      printf("projectpf(time): %lf/%lf/%lf\n",min[5],max[5],sum[5]/np);
      printf("projectpf(iter): %lf/%lf/%lf\n",min[6],max[6],sum[6]/np);
      printf("mg_vcycle      : %lf/%lf/%lf\n",min[7],max[7],sum[7]/np);
      printf("mg_vcycle(down): %lf/%lf/%lf\n",min[8],max[8],sum[8]/np);
      printf("mg_vcycle(up)  : %lf/%lf/%lf\n",min[9],max[9],sum[9]/np);
      printf("Laplacian      : %lf/%lf/%lf\n",min[10],max[10],sum[10]/np);
    }

    GenmapBinSort(h, GENMAP_FIEDLER, &buf0);

    GenmapInt id=GenmapCommRank(local_c);
    int bin;
    if(id<(np+1)/2) bin=0;
    else bin=1;

    GenmapSplitComm(h,&local_c,bin);
    GenmapSetLocalComm(h,local_c);

#if defined(GENMAP_PAUL)
    GenmapBinSort(h,GENMAP_GLOBALID,&buf0);
#endif

    level++;
  }

  crystal_free(&(h->cr));
  buffer_free(&buf0);
}
