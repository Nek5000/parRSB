#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

#define DBG 0

int log2i(sint i){
  sint k=1,l=0;
  while(k<=i) k*=2,l++;
  return l-1;
}

// Following two functions can be combined
void compress_col(struct array *entries){
  GenmapScalar v; sint i,j,je;

  i=j=0; entry *ptr=entries->ptr; while(i<entries->n){
    v=0.0;
    while(j<entries->n && ptr[j].r==ptr[i].r && ptr[j].cn==ptr[i].cn)
      v+=ptr[j].v,j++;
    je=j,j=j-1;
    while(j>=0 && ptr[j].r==ptr[i].r && ptr[j].cn==ptr[i].cn)
      ptr[j].v=v,j--;
    i=j=je;
  }
}

void compress_row(struct array *entries){
  GenmapScalar v; sint i,j,je;

  i=j=0; entry *ptr=entries->ptr; while(i<entries->n){
    v=0.0;
    while(j<entries->n && ptr[j].rn==ptr[i].rn && ptr[j].c==ptr[i].c)
      v+=ptr[j].v,j++;
    je=j,j=j-1;
    while(j>=0 && ptr[j].rn==ptr[i].rn && ptr[j].c==ptr[i].c)
      ptr[j].v=v,j--;
    i=j=je;
  }
}

void mgLevelSetup(parMat M0,mgData data,mgLevel *l_)
{
  uint rn0=M0->rn,nnz0=M0->row_off[rn0];

  struct array entries=null_array;
  array_init(entry,&entries,nnz0);
  entries.n=nnz0;

  uint i,j,nn=0;
  entry *ptr=entries.ptr;
  for(i=0; i<rn0; i++)
    for(j=M0->row_off[i]; j<M0->row_off[i+1]; j++){
      ptr[nn].r=M0->row_start+i;
      ptr[nn].c=M0->col[j];
      ptr[nn].rn=(ptr[nn].r+1)/2;
      ptr[nn].cn=(ptr[nn].c+1)/2; // Let's collapse columns first
      ptr[nn].v=M0->v[j];
      nn++;
    }
  assert(nn==nnz0);

#if DBG
  printf("A: nid=%d nnz=%zu\n",data->c.id,entries.n);
#endif

  slong out[2][1],bf[2][1],in=rn0;
  comm_scan(out,&data->c,gs_long,gs_add,&in,1,bf);
  slong ng=out[1][0];

  ulong ngc=ng/2;
  if(ngc==0) return;

  if(ng>1 && ng%2==1)
    for(j=0; j<M0->row_off[rn0]; j++){
      if(ptr[j].c==ng) ptr[j].cn-=1;
      if(ptr[j].r==ng) ptr[j].rn-=1;
    }

  uint np=data->c.np;
  uint npc=(np/2>0)?(np/2):1;

  if(ngc<npc) npc=ngc;

#if DBG
  printf("id=%u ng=%lld ng_c=%lld np=%d np_c=%d\n",data->c.id,ng,ngc,
      data->c.np,npc);
#endif

  struct crystal cr;
  crystal_init(&cr,&data->c);

  /* coarsen the columns  */
  setOwner(entries.ptr,nnz0,offsetof(entry,cn),offsetof(entry,p),ngc,npc);
  sarray_transfer(entry,&entries,p,1,&cr);

#if DBG
  printf("B: nid=%d nnz=%zu\n",data->c.id,entries.n);
#endif

  buffer buf; buffer_init(&buf,1024);
  if(entries.n){
    sarray_sort_2(entry,entries.ptr,entries.n,r,1,cn,1,&buf);
    compress_col(&entries);
  }

  sarray_transfer(entry,&entries,p,1,&cr);
  assert(entries.n==nnz0);

  /* setup gs ids */
  ptr=entries.ptr;
  slong *ids; GenmapMalloc(rn0,&ids);
  for(i=j=0; i<nnz0; i++)
    if(ptr[i].r==ptr[i].c) ids[j++]=-ptr[i].cn;
  assert(j==rn0);

  /* coarsen the rows */
  setOwner(entries.ptr,nnz0,offsetof(entry,rn),offsetof(entry,p),ngc,npc);
  sarray_transfer(entry,&entries,p,1,&cr);

#if DBG
  printf("C: nid=%d nnz=%zu\n",data->c.id,entries.n);
#endif

  if(entries.n){
    sarray_sort_2(entry,entries.ptr,entries.n,c ,1,rn,1,&buf);
    compress_row(&entries);
    sarray_sort_2(entry,entries.ptr,entries.n,rn,1,cn,1,&buf);
  }

  i=j=nn=0; ptr=entries.ptr; while(i<entries.n){
    while(j<entries.n && ptr[j].rn==ptr[i].rn) j++;
    i=j,nn++;
  }

#if DBG
  printf("D: nid=%d rn=%d\n",data->c.id,nn);
#endif

  /* create the matrix */
  GenmapMalloc(1,l_   ); mgLevel l=*l_ ; l->data=data;
  GenmapMalloc(1,&l->M); parMat M1 =l->M;
  M1->rn=nn; GenmapMalloc(M1->rn+1,&M1->row_off);

  slong cn=nn; comm_scan(out,&data->c,gs_long,gs_add,&cn,1,bf);
  M1->row_start=out[0][0]+1;

  uint nnz1=0;
  i=j=0; ptr=entries.ptr; while(i<entries.n){
    while(j<entries.n && ptr[j].rn==ptr[i].rn && ptr[j].cn==ptr[i].cn)
      j++;
    i=j,nnz1++;
  }

#if DBG
  printf("C: nid=%d rn=%d nnz0=%u nnz1=%d\n",data->c.id,nn,nnz0,nnz1);
#endif

  if(nnz1==0)
    M1->col=NULL,M1->v=NULL;
  else
    GenmapMalloc(nnz1,&M1->col),GenmapMalloc(nnz1,&M1->v);

  uint rn1; GenmapScalar v;
  M1->row_off[0]=i=j=nn=rn1=0; ptr=entries.ptr; while(i<entries.n){
    v=0.0;
    while(j<entries.n && ptr[j].rn==ptr[i].rn && ptr[j].cn==ptr[i].cn)
      j++;
    M1->col[nn]=ptr[i].cn,M1->v[nn]=ptr[i].v,nn++;

    if((j<entries.n && ptr[j].rn!=ptr[i].rn) || j>=entries.n)
      M1->row_off[++rn1]=nn;
    i=j;
  }
  assert(nn==nnz1); //sanity check
  assert(rn1==M1->rn); //sanity check

  GenmapRealloc(rn0+rn1,&ids);
  for(i=nn=0; i<rn1; i++)
    for(j=M1->row_off[i]; j<M1->row_off[i+1]; j++)
      if(M1->row_start+i==M1->col[j]) ids[rn0+nn]=M1->col[j],nn++;
  assert(nn==M1->rn);

  l->J=gs_setup(ids,rn0+M1->rn,&data->c,0,gs_crystal_router,0);

  GenmapFree(ids);
  buffer_free(&buf);
  crystal_free(&cr);
  array_free(&entries);
}

void mgSetup(GenmapComm c,parMat M,mgData *d_){
  GenmapMalloc(1,d_); mgData d=*d_; comm_dup(&d->c,&c->gsc);

  uint np=GenmapCommSize(c); uint rn=M->rn;

  slong out[2][1],bf[2][1],in=rn;
  comm_scan(out,&d->c,gs_long,gs_add,&in,1,bf);
  slong rg=out[1][0];

  d->nLevels=log2i(rg); // Ignoring L_0
#if DBG
  printf("A: nid=%d nLevel=%d\n",d->c.id,d->nLevels);
#endif
  if(d->nLevels==0) return;

  GenmapMalloc(d->nLevels  ,&d->levels   );
  GenmapMalloc(d->nLevels+1,&d->level_off);
  d->level_off[0]=M->rn;

  fflush(stdout);
  uint i; parMat prev=M; parMatPrint(prev,&c->gsc);
  for(i=0;i<d->nLevels;i++){
    mgLevelSetup(prev,d,&d->levels[i]);
    prev=d->levels[i]->M;
    fflush(stdout);
    comm_barrier(&c->gsc);
    parMatPrint(prev,&c->gsc);
    d->level_off[i+1]=d->level_off[i]+prev->rn;
    //TODO: set sigma, nSmooth
  }
  //l0=d->levels[i-1];
  //d->level_off[i]=d->level_off[i-1]+l0->M->row_off[l0->M->rn];

  //GenmapMalloc(d->level_off[i],&d->x  );
  //GenmapMalloc(d->level_off[i],&d->b  );
  //GenmapMalloc(d->level_off[i],&d->buf);
}
