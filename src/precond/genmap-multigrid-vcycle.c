#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

void mg_vcycle(GenmapScalar *u1,GenmapScalar *rhs,mgData d){

  GenmapScalar *s   =d->x;
  GenmapScalar *Gs  =d->y;
  GenmapScalar *r   =d->b;
  GenmapScalar *u   =d->u;

  mgLevel *lvls=d->levels; uint *lvl_off=d->level_off;
  mgLevel l; csr_mat M;

  buffer buf; buffer_init(&buf,1024);

  int nsmooth,nlevels=d->nlevels,lvl;
  GenmapScalar *diag,sigma;
  uint off,n,i,j;

  for(i=0; i<lvl_off[nlevels]; i++)
    s[i]=Gs[i]=r[i]=u[i]=0.0;
  for(i=0; i<lvl_off[1]; i++)
    r[i]=rhs[i];

  GenmapHandle h=d->h;
  comm_barrier(&d->c);
  h->time[8]-=comm_time();

  for(lvl=0; lvl<nlevels-1; lvl++){
    off=lvl_off[lvl]; n=lvl_off[lvl+1]-off;
    l=lvls[lvl]; nsmooth=l->nsmooth; sigma=l->sigma;
    M=l->M; diag=M->diag; assert(n==M->rn);

    //u=sigma*D*rhs
    for(j=0; j<n; j++)
      u[off+j]=sigma*r[off+j]/diag[j];

    // G*u
    comm_barrier(&d->c);
    h->time[10]-=comm_time();
    csr_mat_gather(M,M->gsh,u+off,d->buf,&buf);
    csr_mat_apply(Gs+off,M,d->buf);
    comm_barrier(&d->c);
    h->time[10]+=comm_time();

    // r=rhs-Gu
    for(j=0; j<n; j++)
        r[off+j]=r[off+j]-Gs[off+j];

    for(i=0; i<nsmooth; i++){
      sigma=sigma+0.066666/nsmooth;
      //s=sigma*D*r, u=u+s
      for(j=0; j<n; j++){
        s[off+j] =sigma*r[off+j]/diag[j];
        u[off+j]+=s[off+j];
      }

      //G*s
      comm_barrier(&d->c);
      h->time[10]-=comm_time();
      csr_mat_gather(M,M->gsh,s+off,d->buf,&buf);
      csr_mat_apply(Gs+off,M,d->buf);
      comm_barrier(&d->c);
      h->time[10]+=comm_time();

      //r=r-Gs
      for(j=0; j<n; j++)
        r[off+j]=r[off+j]-Gs[off+j];
    }

    // interpolate to coarser level
    gs(r+off,gs_double,gs_add,1,l->J,&buf);
  }
  comm_barrier(&d->c);
  h->time[8]+=comm_time();

  //coarsest level
  off=lvl_off[nlevels-1]; n=lvl_off[nlevels]-off;

  if(n==1){
    l=lvls[nlevels-1]; M=l->M;
    assert(M->rn==1);
    if(fabs(M->diag[0])>sqrt(GENMAP_TOL))
      u[off]=r[off]/M->diag[0];
    else
      u[off]=0.0;
    r[off]=u[off];
  }

  comm_barrier(&d->c);
  h->time[9]-=comm_time();

  GenmapScalar over=1.33333;
  for(lvl=nlevels-2; lvl>=0; lvl--){
    l=lvls[lvl];
    off=lvl_off[lvl];
    // J*e
    gs(r+off,gs_double,gs_add,0,l->J,&buf);

    //u=u+over*J*e
    n=lvl_off[lvl+1]-off;
    for(j=0; j<n; j++)
      r[off+j]=over*r[off+j]+u[off+j];
  }
  comm_barrier(&d->c);
  h->time[9]+=comm_time();

  // avoid this
  for(i=0; i<lvl_off[1]; i++)
    u1[i]=r[i];

  buffer_free(&buf);
}