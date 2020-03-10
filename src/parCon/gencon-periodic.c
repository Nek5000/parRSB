#include "gencon-impl.h"
#include "exa-memory.h"

#include "math.h"
#include "string.h"
#include "stdio.h"

int faces3D[GC_MAX_FACES][GC_MAX_FACE_VERTICES]={
  {1,5,7,3},{2,4,8,6},{1,2,6,5},{3,7,8,4},{1,3,4,2},{5,6,8,7}
};

int faces2D[GC_MAX_FACES][GC_MAX_FACE_VERTICES]={
  {3,1,0,0},{2,4,0,0},{1,2,0,0},{4,3,0,0},{0,0,0,0},{0,0,0,0}
};

struct minPair_private{
  exaInt proc;
  exaLong orig;
  exaLong min;
};
typedef struct minPair_private* minPair;

int compressPeriodicVertices(exaHandle h,Mesh mesh){
  exaSort(mesh->elements,exaLong_t,\
    offsetof(struct Point_private,globalId),
    exaSortAlgoBinSort,0,exaGetComm(h));

  Point   points=exaArrayGetPointer(mesh->elements);
  exaInt nPoints=exaArrayGetSize(mesh->elements);

  exaInt i,nUnique=0;
  if(nPoints){
    exaLong current=points[0].globalId;
    points[0].globalId=nUnique;
    for(i=1;i<nPoints;i++)
      if(points[i].globalId==current)
        points[i].globalId=nUnique;
      else{
        current=points[i].globalId,++nUnique;
        points[i].globalId=nUnique;
      }
  }

  exaLong out[2][1],buf[2][1],in[1];
  if(nPoints) in[0]=nUnique+1;
  else in[0]=0;

  exaScan(h,out,in,buf,1,exaLong_t,exaAddOp);
  exaLong start=out[0][0];

  for(i=0;i<nPoints;i++) points[i].globalId+=start;

  int size=exaSize(h);
  int rank=exaRank(h);

  exaLong nelgt=mesh->nelgt;
  exaInt nelt=nelgt/size,nrem=nelgt-nelt*size;
  exaLong N=nrem*(nelt+1);

  exaLong eid;
  for(i=0;i<nPoints;i++){
    eid=points[i].elementId;
    if(N==0) points[i].proc=eid/nelt;
    else if(eid+1<=N)
      points[i].proc=ceil((eid+1.0)/(nelt+1.0))-1;
    else points[i].proc=ceil((eid+1.0-N)/nelt)-1+nrem;
  }

  exaComm c=exaGetComm(h);
  exaArrayTransfer(mesh->elements,\
      offsetof(struct Point_private,proc),1,c);
}

exaLong findMinBelowI(exaLong min,exaInt I,exaArray arr){
  minPair ptr=exaArrayGetPointer(arr);

  exaInt i;
  for(i=0;i<I;i++)
    if(ptr[i].orig==min) return ptr[i].min;
  return min;
}

int renumberPeriodicVertices(exaHandle h,Mesh mesh,\
  exaArray matched)
{
  minPair ptr=exaArrayGetPointer(matched);
  exaInt size=exaArrayGetSize(matched);

  exaLong *ids;
  exaMalloc(size,&ids);
  exaInt i;
  for(i=0;i<size;i++) ids[i]=ptr[i].orig;

  exaGS t;
  exaGSSetup(ids,size,exaGetComm(h),0,exaGetDebug(h)>0,&t);

  exaBuffer buf; exaBufferCreate(&buf,1024);
  for(i=0;i<size;i++) ids[i]=ptr[i].min;
  exaGSOp(ids,exaLong_t,exaMinOp,0,t,buf);
  for(i=0;i<size;i++) ptr[i].min=ids[i];
  exaFree(ids);
  exaGSFree(t);

  exaSortArray2(matched,\
      exaLong_t,offsetof(struct minPair_private,orig),\
      exaLong_t,offsetof(struct minPair_private,min ));
  ptr=exaArrayGetPointer(matched);

  exaArray compressed;
  exaArrayInit(&compressed,struct minPair_private,10);

  exaInt sizeMatched=exaArrayGetSize(matched);
  if(sizeMatched) exaArrayAppend(compressed,ptr);
  for(i=1;i<sizeMatched;i++)
    if(ptr[i].orig!=ptr[i-1].orig)
      exaArrayAppend(compressed,&ptr[i]);

  ptr=exaArrayGetPointer(compressed);
  for(i=0;i<exaArrayGetSize(compressed);i++) ptr[i].proc=0;
  exaComm c=exaGetComm(h);
  exaArrayTransfer(compressed,\
      offsetof(struct minPair_private,proc),1,c);

  exaInt rank=exaRank(h);
  if(rank==0){
    exaSortArray2(compressed,\
        exaLong_t,offsetof(struct minPair_private,orig),
        exaLong_t,offsetof(struct minPair_private,min ));
    ptr=exaArrayGetPointer(compressed);
    for(i=0;i<exaArrayGetSize(compressed);i++){
      ptr[i].min=findMinBelowI(ptr[i].min,i,compressed);
    }
  }

  size=exaArrayGetSize(mesh->elements);
  exaInt sizec=0;
  if(rank==0) sizec=exaArrayGetSize(compressed);
  exaMalloc(size+sizec,&ids);

  Point pnt=exaArrayGetPointer(mesh->elements);
  for(i=0;i<size ;i++) ids[i]=pnt[i].globalId;
  for(i=0;i<sizec;i++) ids[size+i]=ptr[i].orig;

  exaGSSetup(ids,size+sizec,exaGetComm(h),0,0,&t);

  for(i=0;i<size ;i++) ids[i]=pnt[i].globalId;
  for(i=0;i<sizec;i++) ids[size+i]=ptr[i].min;
  exaGSOp(ids,exaLong_t,exaMinOp,0,t,buf);
  for(i=0;i<size;i++) pnt[i].globalId=ids[i];

  exaBufferFree(buf);
  exaFree(ids);
  exaArrayFree(compressed);
  exaGSFree(t);
}

int findConnectedPeriodicPairs(exaHandle h,Mesh mesh,\
  BoundaryFace f_,BoundaryFace g_,exaArray matched)
{
  BoundaryFace f,g;
  exaMalloc(1,&f); exaMalloc(1,&g);
  memcpy(f,f_,sizeof(*f)); memcpy(g,g_,sizeof(*g));

  int nvf=mesh->nVertex/2;
  int nDim=mesh->nDim;

  int i,j;
  exaScalar fMax=0.0,gMax=0.0;

  for(i=0;i<nDim;i++){
    exaScalar meanF=0.0,meanG=0.0;

    for(j=0;j<nvf;j++){
      fMax=max(fMax,fabs(f->face.vertex[j].x[i]));
      gMax=max(gMax,fabs(g->face.vertex[j].x[i]));
      meanF+=f->face.vertex[j].x[i];
      meanG+=g->face.vertex[j].x[i];
    }

    for(j=0;j<nvf;j++){
      f->face.vertex[j].x[i]-=(meanF/nvf);
      g->face.vertex[j].x[i]-=(meanG/nvf);
    }
  }

  int shift=0,k;
  exaScalar d2Min=1.e20,d2;
  for(i=0;i<nvf;i++){
    d2=0.0;
    for(j=0;j<nvf;j++){
      k=(j+i)%nvf;
      k=nvf-1-k;
      if(nDim==3)
        d2+=distance3D(f->face.vertex[j],g->face.vertex[k]);
      if(nDim==2)
        d2+=distance2D(f->face.vertex[j],g->face.vertex[k]);
    }
    if(d2<d2Min) {d2Min=d2;shift=i;};
  }
  d2Min=sqrt(d2Min);

  exaScalar fgMax=max(fMax,gMax);
  exaScalar tol=(1e-3)*fgMax;
  if(d2Min>tol){
    fprintf(stderr,"Faces did not match: (d2Min,tol,faceId1"\
        ",faceId2): %lf %lf %lld %lld\n",d2Min,tol,f->faceId,\
        g->faceId);
    exit(1);
  }
  exaDebug(h,"Periodic face match (elementId,faceId): "\
    "(%lld %lld) and (%lld %lld) %lf\n",\
    f->elementId,f->faceId,g->elementId,g->faceId,d2Min);

  struct minPair_private m;
  for(i=0;i<nvf;i++){
    k=(i+shift)%nvf;
    k=nvf-1-k;
    m.min =min(f->face.vertex[i].globalId,\
        g->face.vertex[k].globalId);
    m.orig=max(f->face.vertex[i].globalId,\
        g->face.vertex[k].globalId);
    exaArrayAppend(matched,&m);
  }

  exaFree(f); exaFree(g);
}

int findConnectedPeriodicFaces(exaHandle h,Mesh mesh,\
    exaArray matched)
{
  exaInt bSize=exaArrayGetSize(mesh->boundary);
  BoundaryFace ptr=exaArrayGetPointer(mesh->boundary);
  exaInt i,j;

  for(i=0;i<bSize-1;i++)
    for(j=i+1;j<bSize;j++)
      if(ptr[j].bc[0]==ptr[i].elementId &&\
          ptr[j].bc[1]==ptr[i].faceId){
        findConnectedPeriodicPairs(h,mesh,&ptr[i],&ptr[j],\
            matched);
      }
}

int gatherMatchingPeriodicFaces(exaHandle h,Mesh mesh){
  int size=exaSize(h);
  int rank=exaRank(h);
  BoundaryFace bPtr=exaArrayGetPointer(mesh->boundary);
  int nFaces=exaArrayGetSize(mesh->boundary);

  exaLong nelgt=mesh->nelgt;
  exaInt nelt=nelgt/size,nrem=nelgt-nelt*size;
  exaLong N=nrem*(nelt+1);

  exaInt i; exaLong eid;
  for(i=0;i<nFaces;i++){
    eid=bPtr[i].bc[0];
    if(eid<bPtr[i].elementId){
      exaDebug(h,"Send matching (%d,%d) to (%d,%d).\n",\
          bPtr[i].elementId,bPtr[i].faceId,
          bPtr[i].bc[0]    ,bPtr[i].bc[1]);

      if(N==0) bPtr[i].proc=eid/nelt;
      else if(eid<N) bPtr[i].proc=ceil((eid+1.0)/(nelt+1.0))-1;
      else bPtr[i].proc=ceil((eid+1.0-N)/nelt)-1+nrem;
    } else bPtr[i].proc=rank;
  }

  exaComm c=exaGetComm(h);
  exaArrayTransfer(mesh->boundary,\
      offsetof(struct Boundary_private,proc),1,c);
}

int setPeriodicFaceCoordinates(exaHandle h,Mesh mesh){
  /* Need boundary array to be sorted by elementId */
  exaSortArray(mesh->boundary,exaLong_t,\
      offsetof(struct Boundary_private,elementId));

  BoundaryFace bPtr=exaArrayGetPointer(mesh->boundary);
  exaInt bSize     =exaArrayGetSize(mesh->boundary);
  if(bSize==0) return 0;

  /* Need element array to be sorted by sequenceId */
  exaSortArray(mesh->elements,exaLong_t,\
      offsetof(struct Point_private,sequenceId));

  Point ePtr  =exaArrayGetPointer(mesh->elements);
  exaInt eSize=exaArrayGetSize(mesh->elements);
  if(eSize==0) return 0;

  int faces[GC_MAX_FACES][GC_MAX_FACE_VERTICES];
  if(mesh->nDim==3)
    memcpy(faces,faces3D,\
        GC_MAX_FACES*GC_MAX_FACE_VERTICES*sizeof(int));
  else
    memcpy(faces,faces2D,\
        GC_MAX_FACES*GC_MAX_FACE_VERTICES*sizeof(int));

  exaInt i=0,k=0;
  int nv=mesh->nVertex,nvf=mesh->nVertex/2,j;
  while(i<bSize){
    while(k<eSize && ePtr[k].elementId<bPtr[i].elementId)
      k+=nv;
    //copy vertices to boundary face
    if(k<eSize && ePtr[k].elementId==bPtr[i].elementId){
      int faceId=bPtr[i].faceId;
      for(j=0;j<nvf;j++)
        bPtr[i].face.vertex[j]=ePtr[k+faces[faceId][j]-1];

      exaDebug(h,"Periodic BC (element,face):"
        "(%d,%d)\n",bPtr[i].bc[0],bPtr[i].bc[1]);
    }
    i++;
  }
}

int matchPeriodicFaces(exaHandle h,Mesh mesh){
  setPeriodicFaceCoordinates(h,mesh);
  gatherMatchingPeriodicFaces(h,mesh);

  exaArray matched;
  exaArrayInit(&matched,struct minPair_private,10);

  findConnectedPeriodicFaces(h,mesh,matched);
  renumberPeriodicVertices(h,mesh,matched);
  compressPeriodicVertices(h,mesh);

  exaArrayFree(matched);

  exaSortArray(mesh->elements,exaLong_t,\
      offsetof(struct Point_private,sequenceId));
}
