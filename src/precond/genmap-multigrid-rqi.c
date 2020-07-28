#include <math.h>
#include <stdio.h>

#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

//input r should have zero-sum
int rqi(GenmapHandle h,GenmapComm c,GenmapVector z,int iter,int verbose,
  GenmapVector y)
{
  assert(z->size==y->size);

  GenmapScalar lambda,norm,normi;
  GenmapLong nelg=GenmapGetNGlobalElements(h);

  GenmapVector err; GenmapCreateVector(&err,z->size);

  comm_barrier(&c->gsc);
  h->time[4]-=comm_time();
  mgData d; mgSetup(c,c->M,&d);
  comm_barrier(&c->gsc);
  h->time[4]+=comm_time();

  int rank=GenmapCommRank(c);
  int projecti;

  comm_barrier(&c->gsc);
  h->time[5]-=comm_time();
  h->time[6]+=project_pf(h,c,d,z,30,0,y);
  comm_barrier(&c->gsc);
  h->time[5]+=comm_time();

  lambda=GenmapDotVector(y,z);
  GenmapGop(c,&lambda,1,GENMAP_SCALAR,GENMAP_SUM);

  GenmapAxpbyVector(err,y,1.0,z,-lambda);
  norm=GenmapDotVector(err,err);
  GenmapGop(c,&norm,1,GENMAP_SCALAR,GENMAP_SUM);

  double eps=1e-5;
  GenmapScalar rnorm=norm*eps;

  uint i;
  for(i=0; i<iter && norm>rnorm; i++){
    norm=GenmapDotVector(y,y);
    GenmapGop(c,&norm,1,GENMAP_SCALAR,GENMAP_SUM);
    normi=1.0/sqrt(norm);

    GenmapAxpbyVector(z,z,0.0,y,normi);
    GenmapOrthogonalizebyOneVector(h,c,z,nelg);

    comm_barrier(&c->gsc);
    h->time[5]-=comm_time();
    h->time[6]+=project_pf(h,c,d,z,30,0,y);
    comm_barrier(&c->gsc);
    h->time[5]+=comm_time();

    GenmapOrthogonalizebyOneVector(h,c,y,nelg);

    lambda=GenmapDotVector(y,z);
    GenmapGop(c,&lambda,1,GENMAP_SCALAR,GENMAP_SUM);

    GenmapAxpbyVector(err,y,1.0,z,-lambda);
    norm=GenmapDotVector(err,err);
    GenmapGop(c,&norm,1,GENMAP_SCALAR,GENMAP_SUM);
    if(rank==0 && verbose)
      printf("i=%02d lambda=%1.10e\n",i,lambda);
  }

  mgFree(d);
  GenmapDestroyVector(err);

  return i;
}
