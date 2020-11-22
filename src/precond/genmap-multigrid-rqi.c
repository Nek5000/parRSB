#include <math.h>
#include <stdio.h>

#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

// Input z should be orthogonal to 1-vector, have unit norm.
// RQI should not change z.
int rqi(genmap_handle h,GenmapComm c,mgData d,GenmapVector z,int max_iter,
  int verbose,GenmapVector y)
{
  uint lelt=z->size;
  assert(lelt==y->size);

  struct comm *gsc=&c->gsc;

  metric_tic(gsc,PROJECTPF);
  int ppfi=project_pf(h,c,d,z,20,verbose,y);
  metric_toc(gsc,PROJECTPF);
  metric_acc(NPROJECTPF,ppfi);

  GenmapVector err; GenmapCreateVector(&err,lelt);
  GenmapLong nelg=GenmapGetNGlobalElements(h);
  int rank=GenmapCommRank(GenmapGetGlobalComm(h));

  GenmapScalar *Z,*B,*Bt,*rhs,*solution,*buf;
  GenmapMalloc(max_iter*lelt,&Z);
  GenmapMalloc(max_iter*max_iter,&B);
  GenmapMalloc(max_iter*max_iter,&Bt);
  GenmapMalloc(max_iter*max_iter,&buf);
  GenmapMalloc(max_iter,&rhs);
  GenmapMalloc(max_iter,&solution);

  uint i,j,k,l;
  for(i=0; i<20; i++){
    GenmapScalar norm0=GenmapDotVector(y,y);
    GenmapGop(c,&norm0,1,GENMAP_SCALAR,GENMAP_SUM);
    GenmapScalar normi0=1.0/sqrt(norm0);

    GenmapAxpbyVector(z,z,0.0,y,normi0);
    GenmapOrthogonalizebyOneVector(h,c,z,nelg);

#if defined(GENMAP_GRAMMIAN)
    metric_tic(gsc,GRAMMIAN);
    //Z(:,k) = z;
    //if k>1;
    //  Z(:,k)=z-Z(:,1:k-1)*(B\(Z(:,1:k-1)'*z));
    //  Z(:,k)=Z(:,k)/norm(Z(:,k));
    //end;

    if(i>0){
      // rhs = Z[1:k-1,:]*z
      for(j=0; j<i; j++){
        rhs[j]=0.0;
        for(l=0; l<lelt; l++)
          rhs[j]+=Z[j*lelt+l]*z->data[l];
      }

      //TODO: solution = B\rhs
      // transpose B and invert B and transpose again
      for(j=0; j<i; j++){
        for(k=0; k<i; k++){
          Bt[k*i+j]=B[j*i+k];
        }
      }

      matrix_inverse(i,Bt);

      for(j=0; j<i; j++){
        for(k=0; k<i; k++){
          B[k*i+j]=Bt[j*i+k];
        }
      }

      for(j=0; j<i; j++){
        solution[j]=0.0;
        for(k=0; k<i; k++){
          solution[j]+=B[j*i+k]*rhs[k];
        }
      }

      //Z[k,:] = z[:] - Z[j,:]*solution[j]
      GenmapScalar normz=0.0;
      for(l=0; l<lelt; l++){
        Z[i*lelt+l]=0.0;
        for(j=0; j<i; j++)
          Z[i*lelt+l]+=solution[j]*Z[j*lelt+l];
        Z[i*lelt+l]=z->data[l]-Z[i*lelt+l];
        normz+=Z[i*lelt+l]*Z[i*lelt+l];
      }
      GenmapScalar normzi;
      comm_allreduce(gsc,gs_double,gs_add,&normz,1,&normzi);
      normzi=1.0/normz;

      //Z[k,:]= Z[k,:]/||Z[k,:]||
      for(l=0; l<lelt; l++)
        Z[i*lelt+l]/=normzi;
    }else{
      for(j=0; j<lelt; j++)
        Z[i*lelt+j]=z->data[j];
    }

    //M=Z(:,1:k)'*G*Z(:,1:k);
    //B=Z(:,1:k)'*Z(:,1:k);
    for(j=0; j<i; j++){
      for(k=0; k<i; k++){
        B[j*i+k]=0.0;
        for(l=0; l<lelt; l++)
          B[j*i+k]+=Z[j*lelt+l]*Z[l+k*lelt];
      }
    }

    //global reduction of B
    comm_allreduce(gsc,gs_double,gs_add,B,(i+1)*(i+1),buf);
    metric_toc(gsc,GRAMMIAN);
#endif

    metric_tic(gsc,PROJECTPF);
    ppfi=project_pf(h,c,d,z,20,verbose,y);
    metric_toc(gsc,PROJECTPF);
    metric_acc(NPROJECTPF,ppfi);
    GenmapOrthogonalizebyOneVector(h,c,y,nelg);

    GenmapScalar lambda=GenmapDotVector(y,z);
    GenmapGop(c,&lambda,1,GENMAP_SCALAR,GENMAP_SUM);

    GenmapAxpbyVector(err,y,1.0,z,-lambda);
    GenmapScalar norme=GenmapDotVector(err,err);
    GenmapGop(c,&norme,1,GENMAP_SCALAR,GENMAP_SUM);
    norme=sqrt(norme);

    GenmapScalar norm1=GenmapDotVector(y,y);
    GenmapGop(c,&norm1,1,GENMAP_SCALAR,GENMAP_SUM);
    GenmapScalar normi1=1.0/sqrt(norm1);
  }

  GenmapFree(Z);
  GenmapFree(B);
  GenmapFree(Bt);
  GenmapFree(rhs);
  GenmapFree(solution);
  GenmapFree(buf);

  GenmapDestroyVector(err);

  return i;
}
