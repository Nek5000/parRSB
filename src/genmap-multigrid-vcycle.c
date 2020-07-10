#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

void mg_vcycle(GenmapScalar *u,GenmapScalar *rhs,mgData d){

  GenmapScalar *s   =d->x;
  GenmapScalar *Gs  =d->y;
  GenmapScalar *r   =d->b;

  mgLevel *lvls=d->levels; uint *lvl_off=d->level_off;
  mgLevel l; parMat M;

  buffer buf; buffer_init(&buf,1024);

  int nsmooth,nlevels=d->nlevels;
  GenmapScalar *diag,sigma;
  uint lvl,off,n,i,j;
  for(lvl=0; lvl<nlevels-1; lvl++){
    off=lvl_off[lvl]; n=lvl_off[lvl+1]-off;
    l=lvls[lvl]; nsmooth=l->nsmooth; sigma=l->sigma;
    M=l->M; diag=M->diag;

    //u=sigma*D*rhs
    for(j=0; j<n; j++)
        s[off+j]=sigma*rhs[off+j]/diag[j];

    // r=rhs-G*u
    parMatApply(r+off,M,u+off,d->buf);
    for(j=0; j<n; j++)
        r[off+j]=rhs[off+j]-r[off+j];

    for(i=0; i<nsmooth; i++){
      sigma=sigma+0.066666/nsmooth;
      // s=sigma*D*r, u=u+s
      for(j=0; j<n; j++){
        s[off+j] =sigma*rhs[off+j]/diag[j];
        u[off+j]+=s[off+j];
      }

      //r=r-G*s
      parMatApply(Gs+off,M,s+off,d->buf);
      for(j=0; j<n; j++)
        r[off+j]=r[off+j]-Gs[off+j];

    }
    // interpolate to coarser level
    gs(r+lvl_off[lvl],gs_double,gs_add,1,l->J,&buf);
  }

  //coarsest level
  off=lvl_off[nlevels-1]; n=lvl_off[nlevels]-off;

  if(n==1){
    l=lvls[nlevels-1]; M=l->M;
    assert(M->rn==1);
    if(fabs(M->diag[0])>GENMAP_TOL)
      u[off]=r[off]/M->diag[0];
    else
      u[off]=0.0;
    r[off]=u[off];
  }

  for(lvl=nlevels-2; lvl>=0; lvl++){
    l=lvls[lvl];
    gs(r+lvl_off[lvl+1],gs_double,gs_add,0,l->J,&buf);

    n=lvl_off[lvl+1]-off;
    GenmapScalar over=1.33333;
    for(j=0; j<n; j++)
      r[off+j]=over*r[off+j]+u[off+j];
  }

  // avoid this
  for(i=0; i<lvl_off[nlevels]; i++)
    u[i]=r[i];

  buffer_free(&buf);
}