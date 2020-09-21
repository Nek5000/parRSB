#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

#include <genmap-impl.h>

#define MAXMETS 130
#define MAXLVLS  30
#define MAXSIZE (MAXMETS*MAXLVLS)

static double metrics[MAXMETS];
static double *stack;
static uint stack_size,stack_max;

void metric_init(){
  uint i; for(i=0; i<MAXMETS; i++)
    metrics[i]=0.0;
  stack=NULL;
  stack_size=stack_max=0;
}

void metric_finalize(){
  if(stack!=NULL)
    GenmapFree(stack);
}

void metric_acc(metric m,double count){ metrics[m]+=count; }

void metric_tic(struct comm *c,metric m){
  comm_barrier(c);
  metrics[m]-=comm_time();
}

void metric_toc(struct comm *c,metric m){
  comm_barrier(c);
  metrics[m]+=comm_time();
}

void metric_push_level(){
  assert(stack_size<=stack_max && "stack_size > stack_max");

  if(stack_size==stack_max){
    stack_max+=stack_size/2+1;
    GenmapRealloc(stack_max*MAXMETS,&stack);
  }

  uint i; for(i=0; i<MAXMETS; i++){
    stack[stack_size*MAXMETS+i]=metrics[i];
    metrics[i]=0.0;
  }
  stack_size++;
}

void metric_print(struct comm *c){
  double min[MAXSIZE],max[MAXSIZE],sum[MAXSIZE],buf[MAXSIZE];
  uint max_size=stack_size*MAXMETS;
  assert(max_size<=MAXSIZE);

  uint i; for(i=0; i<max_size; i++)
    min[i]=max[i]=sum[i]=stack[i];

  comm_allreduce(c,gs_double,gs_min,min,MAXSIZE,buf);// min
  comm_allreduce(c,gs_double,gs_max,max,MAXSIZE,buf);// max
  comm_allreduce(c,gs_double,gs_add,sum,MAXSIZE,buf);// sum

  for(i=0; i<max_size; i++)
    sum[i]/=c->np;

#define summary(i,m) sum[i*MAXMETS+m],min[i*MAXMETS+m],max[i*MAXMETS+m]

  for(i=0; i<stack_size; i++){
    if(c->id==0){
      printf("level=%02d\n",i);
#if defined(GENMAP_RCB)
      printf("  RCB             : %g/%g/%g\n",summary(i,RCB));
#endif
      printf("  pairwise        : %g/%g/%g\n",summary(i,PAIRWISE      ));
      printf("  crystal_router  : %g/%g/%g\n",summary(i,CRYSTAL       ));
      printf("  allreduce       : %g/%g/%g\n",summary(i,ALLREDUCE     ));
      printf("  nneighbors      : %g/%g/%g\n",summary(i,NNEIGHBORS    ));
      printf("  gs_setup        : %g/%g/%g\n",summary(i,GSSETUP       ));
      printf("  laplacian_setup : %g/%g/%g\n",summary(i,LAPLACIANSETUP));
      printf("  nconn           : %g/%g/%g\n",summary(i,NCONN         ));
      printf("  gs              : %g/%g/%g\n",summary(i,GSOP          ));
      printf("  laplacian       : %g/%g/%g\n",summary(i,LAPLACIAN     ));
#if defined(GENMAP_RQI) || defined(GENMAP_FMG)
      printf("  precon_setup    : %g/%g/%g\n",summary(i,PRECONSETUP   ));
      printf("  precon_ax       : %g/%g/%g\n",summary(i,PRECONAX      ));
      printf("  precon_vcycle   : %g/%g/%g\n",summary(i,PRECONVCYCLE  ));
      printf("  projectpf       : %g/%g/%g\n",summary(i,PROJECTPF     ));
      printf("  nprojectpf      : %g/%g/%g\n",summary(i,NPROJECTPF    ));
#if defined(GENMAP_RQI)
      printf("  RQI             : %g/%g/%g\n",summary(i,RQI           ));
#elif defined(GENMAP_FMG)
      printf("  FMG             : %g/%g/%g\n",summary(i,FMG           ));
#endif
      for(int j=0; j<100; j++)
      printf("  PPF i=%02d        : %g/%g/%g\n",j,summary(i,END+j       ));
#endif
      printf("  fiedler_time    : %g/%g/%g\n",summary(i,FIEDLER       ));
      printf("  fiedler_iter    : %g/%g/%g\n",summary(i,NFIEDLER      ));
      printf("  RSB             : %g/%g/%g\n",summary(i,RSB           ));
    }
  }
}

#undef summary
#undef MAXMETS
#undef MAXLVLS
#undef MAXSIZE
