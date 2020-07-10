#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

#include <math.h>
#include <stdio.h>

int flex_cg(GenmapHandle h,GenmapComm c,mgData d,GenmapVector r,
  int maxIter,GenmapVector x)
{
  assert(x->size==r->size);
  assert(x->size==GenmapGetNLocalElements(h));

  uint lelt=x->size;

  GenmapVector z0,z,dz,w,p,weights;
  GenmapCreateVector(&z ,lelt); GenmapCreateVector(&w      ,lelt);
  GenmapCreateVector(&p ,lelt); GenmapCreateVector(&weights,lelt);
  GenmapCreateVector(&z0,lelt); GenmapCreateVector(&dz     ,lelt);

  uint i;
  for(i=0; i<lelt; i++) x->data[i]=0,z->data[i]=r->data[i];

  GenmapOrthogonalizebyOneVector(h,c,z,GenmapGetNGlobalElements(h));
  for(i=0; i<lelt; i++) p->data[i]=z->data[i];

  GenmapScalar den,alpha,beta,rz0,rz1,rz2;

  rz1=GenmapDotVector(r,z);
  GenmapGop(c,&rz1,1,GENMAP_SCALAR,GENMAP_SUM);

  GenmapInitLaplacian(h,c,weights);

  for(i=0; i<maxIter; i++){
    GenmapLaplacian(h,c,p,weights,w);

    den=GenmapDotVector(p,w);
    GenmapGop(c,&den,1,GENMAP_SCALAR,GENMAP_SUM);

    alpha=rz1/den;

    GenmapAxpbyVector(x,x,1.0,p, alpha);
    GenmapAxpbyVector(r,r,1.0,w,-alpha);

    GenmapCopyVector(z0,z);
    mg_vcycle(z->data,r->data,d);
    //remove mean from z;

    GenmapAxpbyVector(dz,z,1.0,z0,-1.0);

    rz0=rz1;

    rz1=GenmapDotVector(r,z);
    GenmapGop(c,&rz1,1,GENMAP_SCALAR,GENMAP_SUM);

    rz2=GenmapDotVector(r,dz);
    GenmapGop(c,&rz2,1,GENMAP_SCALAR,GENMAP_SUM);
    
    beta=rz2/rz0;

    GenmapAxpbyVector(p,z,1.0,p,beta);
  }

  GenmapDestroyVector(z) ,GenmapDestroyVector(w      );
  GenmapDestroyVector(p) ,GenmapDestroyVector(weights);
  GenmapDestroyVector(z0),GenmapDestroyVector(dz     );

  return i;
}
