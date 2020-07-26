#include <math.h>
#include <stdio.h>

#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

#define MM 51

int project_pf(GenmapHandle h,GenmapComm c,mgData d,GenmapVector ri,
  int maxIter,int verbose,GenmapVector x)
{
  assert(x->size==ri->size);
  assert(x->size==GenmapGetNLocalElements(h));

  uint lelt=x->size;
  GenmapLong nelg=GenmapGetNGlobalElements(h);

  GenmapVector z0,z,dz,w,p,r;
  GenmapCreateVector(&z ,lelt);
  GenmapCreateVector(&w ,lelt);
  GenmapCreateVector(&r,lelt);
  GenmapCreateVector(&p ,lelt);
  GenmapCreateVector(&z0,lelt);
  GenmapCreateVector(&dz,lelt);

  double *P,*W,a;
  assert(maxIter<MM);
  GenmapCalloc(lelt*MM,&P);
  GenmapCalloc(lelt*MM,&W);

#define PREC 1
#define ORTH 1

  int rank=GenmapCommRank(c);

#if PREC
  if(rank==0 && verbose) printf("Using MG Prec.\n");
#endif
#if ORTH
  if(rank==0 && verbose) printf("Using Orthogonalization.\n");
#endif

  uint i;
  for(i=0; i<lelt; i++)
    x->data[i]=0.0,r->data[i]=ri->data[i];

#if PREC
  mg_vcycle(z->data,r->data,d);
#else
  GenmapCopyVector(z,r);
#endif
#if ORTH
  GenmapOrthogonalizebyOneVector(h,c,z,nelg);
#endif

  GenmapScalar den,alpha,beta,rz0,rz1=0,rz2,rr,scale;

  rz1=GenmapDotVector(r,z);
  GenmapGop(c,&rz1,1,GENMAP_SCALAR,GENMAP_SUM);
  if(GenmapCommRank(c)==0 && verbose)
    printf("rz1=%lf\n",rz1);

  GenmapCopyVector(p,z);

  i=0; uint j,k;
  while(i<maxIter && sqrt(rz1)>GENMAP_TOL){
    GenmapLaplacian(h,c,p,w);
    den=GenmapDotVector(p,w);
    GenmapGop(c,&den,1,GENMAP_SCALAR,GENMAP_SUM);
    alpha=rz1/den;

    scale=1.0/sqrt(den);
    for(j=0; j<lelt; j++){
      W[i*lelt+j]=scale*w->data[j];
      P[i*lelt+j]=scale*p->data[j];
    }

    GenmapAxpbyVector(x,x,1.0,p, alpha);
    GenmapAxpbyVector(r,r,1.0,w,-alpha);

    GenmapCopyVector(z0,z);
#if PREC
    mg_vcycle(z->data,r->data,d);
#else
    GenmapCopyVector(z,r);
#endif
#if ORTH
    GenmapOrthogonalizebyOneVector(h,c,z,nelg);
#endif

    GenmapAxpbyVector(dz,z,1.0,z0,-1.0);

    rr=GenmapDotVector(r,r);
    GenmapGop(c,&rr,1,GENMAP_SCALAR,GENMAP_SUM);
    if(rank==0 && verbose)
      printf("i=%d rr=%1.10e\n",i,sqrt(rr));

    rz0=rz1;

    rz1=GenmapDotVector(r,z);
    GenmapGop(c,&rz1,1,GENMAP_SCALAR,GENMAP_SUM);

    rz2=GenmapDotVector(r,dz);
    GenmapGop(c,&rz2,1,GENMAP_SCALAR,GENMAP_SUM);

    beta=rz2/rz0;

    GenmapAxpbyVector(p,z,1.0,p,beta);

    i++;

    for(k=0; k<lelt; k++)
      P[(MM-1)*lelt+k]=0.0;

    for(j=0; j<i; j++){
      a=0.0;
      for(k=0; k<lelt; k++)
        a+=W[j*lelt+k]*p->data[k];
      GenmapGop(c,&a,1,GENMAP_SCALAR,GENMAP_SUM);
      for(k=0; k<lelt; k++)
        P[(MM-1)*lelt+k]+=a*P[j*lelt+k];
    }

    for(k=0; k<lelt; k++)
      p->data[k]-=P[(MM-1)*lelt+k];
  }

  GenmapDestroyVector(z),GenmapDestroyVector(w);
  GenmapDestroyVector(p); GenmapDestroyVector(r);
  GenmapDestroyVector(z0),GenmapDestroyVector(dz);

  GenmapFree(P); GenmapFree(W);

  return i;
}
