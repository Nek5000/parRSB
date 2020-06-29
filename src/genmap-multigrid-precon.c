#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

int log2i(sint i){
  sint k=1,l=0;
  while(k<=i) k*=2,l++;
  return l-1;
}

void compressMat(struct array *entries,size_t off){
  sint i,i0,j,nn;

  i0=i=j=nn=0; entry *ptr=entries->ptr; while(i<entries->n){
    while(j<entries->n && GETLNG(ptr,j,off)==GETLNG(ptr,i,off))
      if(nn%2) GETLNG(ptr,j,off)=GETLNG(ptr,i-1,off),j++;
      else j++;
    i0=i,i=j,nn++;
  }
  if(nn>1 && nn%2) // take care of last row if its odd
    for(i=i0; i<entries->n; i++) GETLNG(ptr,i,off)=GETLNG(ptr,i0-1,off);
}

void mgLevelSetup(mgLevel l0,mgLevel *l1_)
{
  parMat M0=l0->M; mgData data=l0->data;
  GenmapUInt rn=M0->rn,nnz0=M0->rowOffsets[rn];

  GenmapLong out[2][1],bf[2][1];
  comm_scan(out,&data->c,gs_long,gs_add,&rn,1,bf);
  GenmapLong ng=out[1][0];

  sint np=data->c.np;
  if(ng==np) np/=2;

  printf("ng=%lld np=%d id=%d\n",ng,np,data->c.id);

  struct array entries=null_array;
  array_init(entry,&entries,nnz0);
  entries.n=nnz0;

  GenmapUInt i,j,nn=0;
  entry *ptr=entries.ptr;
  for(i=0; i<rn; i++)
    for(j=M0->rowOffsets[i]; j<M0->rowOffsets[i+1]; j++){
      ptr[nn].r=ptr[nn].rn=i+M0->rowStart ;
      ptr[nn].c=ptr[nn].cn=  M0->colIdx[j];
      ptr[nn].v=M0->v[j],ptr[nn].p=M0->owner[j];
      nn++;
    }
  assert(nn==nnz0);
  printf("A: nid=%d n=%zu\n",data->c.id,entries.n);

  struct crystal cr;
  crystal_init(&cr,&data->c);

  /* coarsen the columns  */
  setOwner(entries.ptr,nnz0,offsetof(entry,c),offsetof(entry,p),ng,np);
  sarray_transfer(entry,&entries,p,1,&cr);

  printf("B: nid=%d n=%zu\n",data->c.id,entries.n);

  buffer buf; buffer_init(&buf,1024);
  if(entries.n){
    sarray_sort_2(entry,entries.ptr,entries.n,c,1,r,1,&buf);
    compressMat(&entries,offsetof(entry,cn));
  }

  sarray_transfer(entry,&entries,p,1,&cr);
  assert(entries.n==nnz0);

  /* setup gs ids */
  ptr=entries.ptr;
  GenmapLong *ids; GenmapMalloc(nnz0,&ids);
  for(i=0; i<nnz0; i++){
    ids[i]=-ptr[i].cn;
    printf("nid=%d nnz=%u cn=%lld\n",data->c.id,nnz0,ptr[i].cn);
  }

  /* coarsen the rows */
  setOwner(entries.ptr,nnz0,offsetof(entry,r),offsetof(entry,p),ng,np);
  sarray_transfer(entry,&entries,p,1,&cr);

  if(entries.n){
    sarray_sort_2(entry,entries.ptr,entries.n,r,1,c,1,&buf);
    compressMat(&entries,offsetof(entry,rn));
  }
  ptr=entries.ptr;

#define TMPDBG 0
#if TMPDBG
  sarray_transfer(entry,&entries,p,0,&cr);
  ptr=entries.ptr;
  for(i=0; i<nnz0; i++)
    printf("nid=%d nnz=%d rn=%lld\n",data->c.id,nnz0,ptr[i].cn);
#endif

  /* create the matrix */
  GenmapMalloc(1,l1_   ); mgLevel l1=*l1_; l1->data=data;
  GenmapMalloc(1,&l1->M); parMat M1 =l1->M; M1->rank=data->c.id;

  sarray_sort_2(entry,ptr,entries.n,rn,1,cn,1,&buf);
  i=j=nn=0; ptr=entries.ptr; while(i<entries.n){
    while(j<entries.n && ptr[j].rn==ptr[i].rn) j++;
    i=j,nn++;
  }
  M1->rn=nn;

  GenmapLong cn=nn; comm_scan(out,&data->c,gs_long,gs_add,&cn,1,bf);
  M1->cn=out[1][0];
  M1->rowStart=out[0][0]+1;

  GenmapUInt nnz1=0;
  i=j=0; ptr=entries.ptr; while(i<entries.n){
    while(j<entries.n && ptr[j].rn==ptr[i].rn && ptr[j].cn==ptr[i].cn)
      j++;
    i=j,nnz1++;
  }
  printf("C: nid=%d nrows=%d nnz0=%u nnz1=%d\n",data->c.id,nn,nnz0,nnz1);

  GenmapMalloc(M1->rn+1,&M1->rowOffsets);
  if(nnz1==0)
    M1->colIdx=NULL,M1->v=NULL,M1->owner=NULL;
  else
    GenmapMalloc(nnz1,&M1->colIdx),GenmapMalloc(nnz1,&M1->v),\
      GenmapMalloc(nnz1,&M1->owner);

  GenmapScalar v;
  sint nr=0; M1->rowOffsets[nr]=0;
  i=j=nn=0; ptr=entries.ptr; while(i<entries.n){
    v=0.0;
    while(j<entries.n && ptr[j].rn==ptr[i].rn && ptr[j].cn==ptr[i].cn)
      v+=ptr[j].v,j++;
    M1->colIdx[nn]=ptr[i].cn,M1->v[nn]=v,nn++;
    //TODO: set owner -- M1->owner[nn]=ptr[i].p,

    if((j<entries.n && ptr[j].rn!=ptr[i].rn) || j>=entries.n)
      M1->rowOffsets[nr++]=nn;
    i=j;
  }
  assert(nn==nnz1); //sanity check
  assert(nr==M1->rn); //sanity check

  GenmapRealloc(nnz0+nnz1,&ids);
  for(i=nnz0; i<nnz0+nnz1; i++)
    if(M1->owner[i]!=M1->rank)
      ids[i]=-M1->colIdx[i-nnz0];
    else
      ids[i]=M1->colIdx[i-nnz0];

  M1->gsh=gs_setup(ids,nnz0+nnz1,&data->c,0,gs_crystal_router,0);

  GenmapFree(ids);
  buffer_free(&buf);
  crystal_free(&cr);
  array_free(&entries);
}

void mgSetup(GenmapComm c,parMat M,mgData *d_){
  GenmapMalloc(1,d_); mgData d=*d_; comm_dup(&d->c,&c->gsc);

  sint np=GenmapCommSize(c);
  sint rn=M->rn;
  d->nLevels=log2i(rn)+log2i(np)+1; // +1 for L_0

  printf("A: nLevel=%d rn=%d np=%d\n",d->nLevels,log2i(rn),log2i(np));

  GenmapGop(c,(void*)&d->nLevels,1,GENMAP_INT,GENMAP_MAX);
  if(d->nLevels==0) return;

  printf("B: nLevel=%d\n",d->nLevels);

  GenmapMalloc(d->nLevels,&d->levels);
  GenmapMalloc(1,&d->levels[0]); d->levels[0]->M=M; d->levels[0]->data=d;
  //TODO: set sigma, nSmooth

  GenmapMalloc(d->nLevels,&d->levelOffsets);
  d->levelOffsets[0]=0;

  sint i; mgLevel l0;
  for(i=1;i<d->nLevels;i++){
    l0=d->levels[i-1];
    d->levelOffsets[i]=d->levelOffsets[i-1]+l0->M->rowOffsets[l0->M->rn];
    mgLevelSetup(l0,&d->levels[i]);
    //TODO: set sigma, nSmooth
  }
  //l0=d->levels[i-1];
  //d->levelOffsets[i]=d->levelOffsets[i-1]+l0->M->rowOffsets[l0->M->rn];

  //GenmapMalloc(d->levelOffsets[i],&d->x  );
  //GenmapMalloc(d->levelOffsets[i],&d->b  );
  //GenmapMalloc(d->levelOffsets[i],&d->buf);
}
