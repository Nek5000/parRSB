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

#define PREC 0
#define LAPO 0
#define ORTH 0

  uint i;
  for(i=0; i<lelt; i++)
    x->data[i]=0.0;
#if PREC
  printf("Using flex\n");
  mg_vcycle(z->data,r->data,d);
#if ORTH
  GenmapOrthogonalizebyOneVector(h,c,z,nelg);
#endif
#else
  GenmapCopyVector(z,r);
#endif

  GenmapScalar den,alpha,beta,rz0,rz1,rz2;

  rz1=GenmapDotVector(r,z);
  GenmapGop(c,&rz1,1,GENMAP_SCALAR,GENMAP_SUM);

  GenmapCopyVector(p,z);

  i=0;
  while(i<maxIter && sqrt(rz1)>GENMAP_TOL){
#if LAPO
    printf("Using original Laplacian\n");
    GenmapLaplacian(h,c,p,weights,w);
#else
    printf("Using weighted Laplacian\n");
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
#if ORTH
    GenmapOrthogonalizebyOneVector(h,c,z,nelg);
#endif
#else
    GenmapCopyVector(z,r);
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
