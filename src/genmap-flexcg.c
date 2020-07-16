#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

#include <math.h>
#include <stdio.h>

/* Orthogonalize by 1-vector (vector of all 1's) */
int ortho_one_vector(GenmapHandle h,GenmapComm c,GenmapVector q1,
  GenmapLong n)
{
  GenmapInt i;
  GenmapScalar sum = 0.0;
  for(i = 0;  i < q1->size; i++) {
    sum += q1->data[i];
  }

  GenmapGop(c, &sum, 1, GENMAP_SCALAR, GENMAP_SUM);
  sum /= n;
  for(i = 0;  i < q1->size; i++) {
    q1->data[i] -= sum;
  }

  return 0;
}

int flex_cg(GenmapHandle h,GenmapComm c,mgData d,GenmapVector r,
  GenmapVector weights,int maxIter,GenmapVector x)
{
  assert(x->size==r->size);
  assert(x->size==GenmapGetNLocalElements(h));

  uint lelt=x->size;
  GenmapLong nelg=GenmapGetNGlobalElements(h);

  GenmapVector z0,z,dz,w,p;
  GenmapCreateVector(&z ,lelt);
  GenmapCreateVector(&w ,lelt);
  GenmapCreateVector(&p ,lelt);
  GenmapCreateVector(&z0,lelt);
  GenmapCreateVector(&dz,lelt);

#define PREC 1
#define ORTH 1
#define LAPO 1

  int rank=GenmapCommRank(c);
#if LAPO
  if(rank==0) printf("Using original Laplacian.\n");
#else
  if(rank==0) printf("Using weighted Laplacian.\n");
#endif

#if PREC
  if(rank==0) printf("Using MG Prec.\n");
#endif

  uint i;
  for(i=0; i<lelt; i++)
    x->data[i]=0.0;
#if PREC
  mg_vcycle(z->data,r->data,d);
#else
  GenmapCopyVector(z,r);
#endif
#if ORTH
  GenmapOrthogonalizebyOneVector(h,c,z,nelg);
#endif

  GenmapScalar den,alpha,beta,rz0,rz1,rz2;

  rz1=GenmapDotVector(r,z);
  GenmapGop(c,&rz1,1,GENMAP_SCALAR,GENMAP_SUM);

  GenmapCopyVector(p,z);

  i=0;
  while(i<maxIter && sqrt(rz1)>1e-10){
#if LAPO
    GenmapLaplacian(h,c,p,weights,w);
#else
    GenmapLaplacianWeighted(h,c,p,weights,w);
#endif

    den=GenmapDotVector(p,w);
    GenmapGop(c,&den,1,GENMAP_SCALAR,GENMAP_SUM);

    alpha=rz1/den;

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

    rz0=rz1;

    rz1=GenmapDotVector(r,z);
    GenmapGop(c,&rz1,1,GENMAP_SCALAR,GENMAP_SUM);

    GenmapAxpbyVector(dz,z,1.0,z0,-1.0);
    rz2=GenmapDotVector(r,dz);
    GenmapGop(c,&rz2,1,GENMAP_SCALAR,GENMAP_SUM);
    beta=rz2/rz0;

    GenmapAxpbyVector(p,z,1.0,p,beta);
    i++;
  }

  GenmapDestroyVector(z),GenmapDestroyVector(w);
  GenmapDestroyVector(p);
  GenmapDestroyVector(z0),GenmapDestroyVector(dz);

  return i;
}
