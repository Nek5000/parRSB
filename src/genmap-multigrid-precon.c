#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

int log2i(GenmapInt i){
  GenmapInt k=1,l=0;
  while(k<=i) k*=2,l++;
  return l-1;
}

void combineRows(GenmapInt in,GenmapLong *icol,GenmapScalar *ival,
    GenmapInt *on,GenmapLong *ocol,GenmapScalar *oval)
{
  typedef struct{GenmapLong c; GenmapScalar v;} entry;

  struct array entries=null_array;
  array_init(entry,&entries,in);
  entries.n=in;

  GenmapInt i; entry *ptr=entries.ptr;
  for(i=0; i<in; i++) ptr[i].c=icol[i],ptr[i].v=ival[i];

  buffer buf; buffer_init(&buf,1024);
  sarray_sort(entry,ptr,in,c,1,&buf);
  buffer_free(&buf);

  GenmapInt j=0,nn=0; GenmapScalar v;
  i=0;
  while(i<in){
    v=ptr[i].v,j=i+1;
    while(j<in && ptr[i].c==ptr[j].c) v+=ptr[j].v,j++;
    ocol[nn]=ptr[i].c,oval[nn]=v;
    i=j,nn++;
  }
  *on=nn;
}

void mgLevelSetup(mgLevel l0,mgLevel *l1_)
{
  //FIXME: Assumes the 'reverse-Nek' element distribution
  parMat M0=l0->M;

  GenmapMalloc(1,l1_   ); mgLevel l1=*l1_;
  GenmapMalloc(1,&l1->M); parMat  M1=l1->M;

  M1->rn=M0->rn/2;
  M1->rank=l0->c.id;

  GenmapLong in[1],out[2][1],buf[2][1];
  in[0]=M1->rn;
  comm_scan(out,&l0->c,gs_long,gs_add,in,1,buf);

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
    for(i=0; i<end; i++){ // combine row 2*i and 2*i+1 in M0
      in=M0->rowOffsets[2*i+2]-M0->rowOffsets[2*i];
      M1->rowOffsets[M1->rn]+=in;
      GenmapRealloc(M1->rowOffsets[M1->rn],&M1->colIdx);
      GenmapRealloc(M1->rowOffsets[M1->rn],&M1->v);
      GenmapRealloc(M1->rowOffsets[M1->rn],&M1->owner);
      combineRows(
        in  ,M0->colIdx+M0->rowOffsets[2*i],M0->v+M0->rowOffsets[2*i],
        &out,M1->colIdx+M1->rowOffsets[i]  ,M0->v+M0->rowOffsets[i]
      );
      M1->rowOffsets[M1->rn]=M1->rowOffsets[i+1]=M1->rowOffsets[i]+out;
    }

    if(odd){ // do the residual
      in=M0->rowOffsets[M0->rn]-M0->rowOffsets[M0->rn-3];
      M1->rowOffsets[M1->rn]+=in;
      GenmapRealloc(M1->rowOffsets[M1->rn],&M1->colIdx);
      GenmapRealloc(M1->rowOffsets[M1->rn],&M1->v);
      GenmapRealloc(M1->rowOffsets[M1->rn],&M1->owner);
      combineRows(in,M0->colIdx+M0->rowOffsets[M0->rn-3],
          M0->v+M0->rowOffsets[M0->rn-3],
        &out,M1->colIdx+M1->rowOffsets[end+1],M0->v+M0->rowOffsets[end+1]
      );
      assert(end+1==M1->rn);
      M1->rowOffsets[M1->rn]=M1->rowOffsets[i]+out;
    }
  }
}

void mgSetup(GenmapComm c,parMat M,mgData *d_){
  GenmapMalloc(1,d_); mgData d=*d_;

  GenmapInt np=GenmapCommSize(c);
  GenmapInt rn=M->rn;
  d->nLevels=log2i(rn)+log2i(np);

  GenmapGop(c,(void*)&d->nLevels,1,GENMAP_INT,GENMAP_MAX);

  if(d->nLevels==0) return;

  GenmapMalloc(d->nLevels,&d->levels);
  d->levels[0]->M=M; comm_dup(&d->levels[0]->c,&c->gsc);
  //TODO: set sigma, nSmooth

  GenmapInt i;
  for(i=1;i<d->nLevels;i++){
    mgLevelSetup(d->levels[i-1],&d->levels[i]);
    // set level offsets
  }
  // allocate d->x and d->b
}
