#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

int log2i(GenmapInt i){
  GenmapInt k=1,l=0;
  while(k<=i) k*=2,l++;
  return l-1;
}

typedef struct{GenmapLong r,c,rn,cn; GenmapInt p; GenmapScalar v;} entry;

#define GETLNG(ptr,i,offset)\
  (*((GenmapLong*)((char*)(ptr)+offset+(i)*sizeof(entry))))

void compressMat(struct array *entries,size_t off){
  GenmapInt i,i0,j,nn;

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

  GenmapInt np=data->c.np;
  if(ng==np) np/=2;

  struct array entries=null_array;
  array_init(entry,&entries,nnz0);
  entries.n=nnz0;

  GenmapUInt i,j,nn=0;
  entry *ptr=entries.ptr;
  for(i=0; i<rn; i++)
    for(j=M0->rowOffsets[i]; j<M0->rowOffsets[i+1]; j++){
      ptr[nn].r=ptr[nn].rn=i+M0->rowStart;
      ptr[nn].c=ptr[nn].cn=M0->colIdx[j];
      ptr[nn].v=M0->v[j],ptr[nn].p=M0->owner[j];
      nn++;
    }
  assert(nn==nnz0);

  //allocate space for gs ids
  GenmapScalar *ids; GenmapMalloc(nnz0,&ids);

  struct crystal cr;
  crystal_init(&cr,&data->c);

  /* coarsen the columns  */
  setOwner(entries.ptr,nnz0,offsetof(entry,c),offsetof(entry,p),ng,np);
  sarray_transfer(entry,&entries,p,1,&cr);

  buffer buf; buffer_init(&buf,1024);
  sarray_sort_2(entry,ptr,entries.n,c,1,r,1,&buf);
  compressMat(&entries,offsetof(entry,cn));

  sarray_transfer(entry,&entries,p,1,&cr);

  /* setup gsids */
  for(i=0; i<nnz0; i++)
    ids[i]=-ptr[i].cn;

  /* coarsen the rows */
  setOwner(entries.ptr,nnz0,offsetof(entry,r),offsetof(entry,p),ng,np);
  sarray_transfer(entry,&entries,p,1,&cr);

  sarray_sort_2(entry,ptr,entries.n,r,1,c,1,&buf);
  compressMat(&entries,offsetof(entry,rn));

  sarray_transfer(entry,&entries,p,0,&cr);

  /* create the matrix */
  GenmapMalloc(1,l1_   ); mgLevel l1=*l1_; l1->data=data;
  GenmapMalloc(1,&l1->M); parMat M1 =l1->M;

  sarray_sort_2(entry,ptr,entries.n,r,1,c,1,&buf);
  i=j=nn=0; ptr=entries.ptr; while(i<entries.n){
    while(j<entries.n && ptr[j].rn==ptr[i].rn) j++;
    i=j,nn++;
  }

  GenmapLong cn=nn;
  comm_scan(out,&data->c,gs_long,gs_add,&nn,1,bf);
  M1->rn=nn,M1->cn=out[1][0];
  M1->rank=data->c.id;
  M1->rowStart=out[0][0]+1;

  GenmapMalloc(M1->rn+1,&M1->rowOffsets);

  GenmapUInt nnz1=0;
  i=j=0; ptr=entries.ptr; while(i<entries.n){
    while(j<entries.n && ptr[j].rn==ptr[i].rn && ptr[j].cn==ptr[i].cn)
      j++;
    i=j,nnz1++;
  }

  if(nnz1==0)
    M1->colIdx=NULL,M1->v=NULL,M1->owner=NULL;
  else
    GenmapMalloc(nnz1,&M1->colIdx),GenmapMalloc(nnz1,&M1->v),\
      GenmapMalloc(nnz1,&M1->owner);

  GenmapScalar v;
  GenmapInt nr;
  i=j=nn=nr=0; ptr=entries.ptr; while(i<entries.n){
    v=0.0,M1->rowOffsets[nr++]=i;
    while(j<entries.n && ptr[j].rn==ptr[i].rn && ptr[j].cn==ptr[i].cn)
      v+=ptr[j].v,j++;
    M1->colIdx[nn]=ptr[i].cn,M1->v[nn]=v,M1->owner[nn]=ptr[i].p,nn++;

    if(j<entries.n && ptr[j].rn!=ptr[i].rn) M1->rowOffsets[nr++]=j;
    i=j;
  }
  assert(nn==nnz1); //sanity check
  assert(nr==M1->rn); //sanity check
  M1->rowOffsets[nr]=nnz1;

  GenmapRealloc(nnz0+nnz1,&ids);
  for(i=nnz0; i<nnz0+nnz1; i++)
    ids[i]=M1->colIdx[i-nnz0];

  //M1->gsh=gs_setup();

  GenmapFree(ids);
  buffer_free(&buf);
  crystal_free(&cr);
  array_free(&entries);
}

void mgSetup(GenmapComm c,parMat M,mgData *d_){
  GenmapMalloc(1,d_); mgData d=*d_;

  GenmapInt np=GenmapCommSize(c);
  GenmapInt rn=M->rn;
  d->nLevels=log2i(rn)+log2i(np);

  GenmapGop(c,(void*)&d->nLevels,1,GENMAP_INT,GENMAP_MAX);

  if(d->nLevels==0) return;

  GenmapMalloc(d->nLevels,&d->levels);
  d->levels[0]->M=M; comm_dup(&d->c,&c->gsc);
  //TODO: set sigma, nSmooth

  GenmapMalloc(d->nLevels,&d->levelOffsets);
  d->levelOffsets[0]=0;

  GenmapInt i; mgLevel l0;
  for(i=1;i<d->nLevels;i++){
    l0=d->levels[i-1];
    d->levelOffsets[i]=d->levelOffsets[i-1]+l0->M->rowOffsets[l0->M->rn];
    mgLevelSetup(l0,&d->levels[i]);
    //TODO: set sigma, nSmooth
  }
  l0=d->levels[i-1];
  d->levelOffsets[i]=d->levelOffsets[i-1]+l0->M->rowOffsets[l0->M->rn];

  GenmapMalloc(d->levelOffsets[i],&d->x  );
  GenmapMalloc(d->levelOffsets[i],&d->b  );
  GenmapMalloc(d->levelOffsets[i],&d->buf);
}

#undef GETPTR

#if 0
void mgLevelSetup(mgLevel l0,mgLevel *l1_)
{
  //FIXME: Assumes the 'reverse-Nek' element distribution
  parMat M0=l0->M;

  GenmapMalloc(1,l1_   ); mgLevel l1=*l1_;
  GenmapMalloc(1,&l1->M); parMat  M1=l1->M;

  M1->rn=M0->rn/2;
  M1->rank=l0->data->c.id;

  GenmapLong in[1],out[2][1],buf[2][1];
  in[0]=M1->rn;
  comm_scan(out,&l0->data->c,gs_long,gs_add,in,1,buf);

  M1->cn      =out[1][0];
  M1->rowStart=out[0][0]+1;

  GenmapMalloc(M1->rn+1,&M1->rowOffsets);
  M1->rowOffsets[0]=M1->rowOffsets[M1->rn]=0;

  if(M1->rn==0){ // or equivalently M0->rn==1
    //regular setup is over, need to setup new comm.
  }else{
    GenmapInt end=M1->rn;
    int odd=M0->rn%2;
    if(odd) end=end-1;

    GenmapInt i,in,out;
    GenmapUInt ooff,ioff;
    for(i=0; i<end; i++){ // combine row 2*i and 2*i+1 in M0
      in=M0->rowOffsets[2*i+2]-M0->rowOffsets[2*i];
      M1->rowOffsets[M1->rn]+=in;
      GenmapRealloc(M1->rowOffsets[M1->rn],&M1->colIdx);
      GenmapRealloc(M1->rowOffsets[M1->rn],&M1->v);
      GenmapRealloc(M1->rowOffsets[M1->rn],&M1->owner);

      ioff=M0->rowOffsets[2*i];
      ooff=M1->rowOffsets[i  ];
      combineRows(
        in  ,M0->colIdx+ioff,M0->v+ioff,M0->owner+ioff,
        &out,M1->colIdx+ooff,M0->v+ooff,M0->owner+ooff
      );
      M1->rowOffsets[M1->rn]=M1->rowOffsets[i+1]=M1->rowOffsets[i]+out;
    }

    if(odd){ // do the residual
      in=M0->rowOffsets[M0->rn]-M0->rowOffsets[M0->rn-3];
      M1->rowOffsets[M1->rn]+=in;
      GenmapRealloc(M1->rowOffsets[M1->rn],&M1->colIdx);
      GenmapRealloc(M1->rowOffsets[M1->rn],&M1->v);
      GenmapRealloc(M1->rowOffsets[M1->rn],&M1->owner);

      ioff=M0->rowOffsets[M0->rn-3];
      ooff=M1->rowOffsets[end];
      combineRows(
        in  ,M0->colIdx+ioff,M0->v+ioff,M0->owner+ioff,
        &out,M1->colIdx+ooff,M0->v+ooff,M1->owner+ooff
      );
      M1->rowOffsets[M1->rn]=M1->rowOffsets[i]+out;
      assert(end+1==M1->rn);
    }
  }

  //compress over columns
}

void combineRows(
    GenmapInt in ,GenmapLong *icol,GenmapScalar *ival,GenmapInt *iowners,
    GenmapInt *on,GenmapLong *ocol,GenmapScalar *oval,GenmapInt *oowners)
{
  typedef struct{GenmapLong c; GenmapInt owner; GenmapScalar v;} entry;

  struct array entries=null_array;
  array_init(entry,&entries,in);
  entries.n=in;

  GenmapInt i; entry *ptr=entries.ptr;
  for(i=0; i<in; i++)
    ptr[i].c=icol[i],ptr[i].v=ival[i],ptr[i].owner=iowners[i];

  buffer buf; buffer_init(&buf,1024);
  sarray_sort(entry,ptr,in,c,1,&buf);
  buffer_free(&buf);

  GenmapInt j=0,nn=0; GenmapScalar v;
  i=0;
  while(i<in){
    v=ptr[i].v,j=i+1;
    while(j<in && ptr[i].c==ptr[j].c) v+=ptr[j].v,j++;
    ocol[nn]=ptr[i].c,oval[nn]=v,oowners[nn]=ptr[i].owner;
    i=j,nn++;
  }
  *on=nn;

  array_free(&entries);
}

void combineCols(
    GenmapInt in ,GenmapLong *icol,GenmapScalar *ival,GenmapInt *iowners,
    GenmapInt *on,GenmapLong *ocol,GenmapScalar *oval,GenmapInt *oowners)
{
  GenmapInt j=0,nn=0,p=0; GenmapScalar v;
  i=0;
  while(i<in){
    p=iowners[i];
    while(j<in && iowners[j]==iowners[i]){ // same proc
      v+=ptr[j].v,j++;
    }
    ocol[nn]=ptr[i].c,oval[nn]=v,oowners[nn]=ptr[i].owner;
    i=j,nn++;
  }
  *on=nn;
}
#endif
