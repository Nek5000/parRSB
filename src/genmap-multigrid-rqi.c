#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

#include <math.h>
#include <stdio.h>

//input r should have zero-sum
int rqi(GenmapHandle h,GenmapComm c,GenmapVector z,int iter,int verbose,
    GenmapVector y)
{
  assert(z->size==y->size);

  GenmapScalar lambda,norm,normi;
  GenmapLong nelg=GenmapGetNGlobalElements(h);

  GenmapVector err; GenmapCreateVector(&err,z->size);
  mgData d; mgSetup(c,c->M,&d);

  int rank=GenmapCommRank(c);

  project_pf(h,c,d,z,30,0,y);

  norm=GenmapDotVector(y,y);
  GenmapGop(c,&norm,1,GENMAP_SCALAR,GENMAP_SUM);
  if(rank==0 && verbose)
    printf("norm(y): %lf\n",sqrt(norm));

  lambda=GenmapDotVector(y,z);
  GenmapGop(c,&lambda,1,GENMAP_SCALAR,GENMAP_SUM);
  if(rank==0 && verbose)
    printf("lambda: %lf\n",lambda);

  uint i;
  for(i=0; i<iter; i++){
    norm=GenmapDotVector(y,y);
    GenmapGop(c,&norm,1,GENMAP_SCALAR,GENMAP_SUM);
    normi=1.0/sqrt(norm);

    GenmapAxpbyVector(z,z,0.0,y,normi);
    GenmapOrthogonalizebyOneVector(h,c,z,nelg);

    project_pf(h,c,d,z,30,0,y);
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
