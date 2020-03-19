#include <math.h>
#include <string.h>

#include "gencon-impl.h"
#include "exa-memory.h"

#define readT(coords,buf,T,nVertex) do{\
  memcpy(coords,buf,sizeof(T)*nVertex);\
} while(0)

#define writeInt(dest,val) do{\
  memcpy(dest,&(val),sizeof(int));\
} while(0)

int PRE_TO_SYM_VERTEX[GC_MAX_VERTICES]={0,1,3,2,4,5,7,6};
int PRE_TO_SYM_FACE[GC_MAX_FACES]={2,1,3,0,4,5};

int transferBoundaryFaces(exaHandle h,Mesh mesh){
  int size=exaSize(h);
  BoundaryFace ptr=exaArrayGetPointer(mesh->boundary);
  int nFaces=exaArrayGetSize(mesh->boundary);

  exaLong nelgt=mesh->nelgt;
  exaInt nelt=nelgt/size,nrem=nelgt-nelt*size;
  exaLong N=nrem*(nelt+1);

  exaInt i; exaLong eid;
  for(i=0;i<nFaces;i++,ptr++){
    eid=ptr->elementId;
    if(N==0) ptr->proc=eid/nelt;
    else if(eid+1<=N) ptr->proc=ceil((eid+1.0)/(nelt+1.0))-1;
    else ptr->proc=ceil((eid+1.0-N)/nelt)-1+nrem;
  }

  exaComm c=exaGetComm(h);
  exaArrayTransfer(mesh->boundary,
    offsetof(struct Boundary_private,proc),1,c);
}

int readRe2Header(exaHandle h,Mesh mesh,MPI_File file){
  char version[BUFSIZ];
  int nelgt,nDim,nelgv,err;
  MPI_Status st;

  exaInt rank=exaRank(h);
  exaInt size=exaSize(h);
  MPI_Comm comm=exaGetMPIComm(h);

  char *buf=(char *)calloc(GC_RE2_HEADER_LEN+1,sizeof(char));
  if(rank==0){
    err=MPI_File_read(file,buf,GC_RE2_HEADER_LEN,MPI_BYTE,&st);
    if(err) return 1;
  }
  MPI_Bcast(buf,GC_RE2_HEADER_LEN,MPI_BYTE,0,comm);
  exaDebug(h,"header: %s\n",buf);
  sscanf(buf,"%5s %d %d %d",version,&nelgt,&nDim,&nelgv);

  int nVertex=(nDim==2)?4:8;
  int nelt=nelgt/size,nrem=nelgt-nelt*size;
  nelt+=(rank<nrem ? 1: 0);

  // Initialize the mesh structure
  mesh->nDim=nDim;
  mesh->nVertex=nVertex;
  mesh->nNeighbors=nDim;
  mesh->nelgt=nelgt;
  mesh->nelgv=nelgv;
  mesh->nelt=nelt;

  if(rank==0){
    exaDebug(h,"ndim,nvertex,nneighbors,nelgt,nelt:%d %d %d %d "
      "%d\n",mesh->nDim,mesh->nVertex,mesh->nNeighbors,\
      mesh->nelgt,mesh->nelt);
  }

  free(buf);

  return 0;
}

int readRe2Coordinates(exaHandle h,Mesh mesh,MPI_File file){
  exaInt rank=exaRank(h);
  exaInt size=exaSize(h);
  MPI_Comm comm=exaGetMPIComm(h);

  int nelt=mesh->nelt;
  int nelgt=mesh->nelgt;
  int nDim=mesh->nDim;
  int nVertex=mesh->nVertex;

  exaLong out[2][1],buff[2][1],in[1];
  in[0]=nelt;
  exaScan(h,out,in,buff,1,exaLong_t,exaAddOp);
  exaLong start=out[0][0];

  int elemDataSize=nVertex*nDim*sizeof(double)+sizeof(double);
  int headerSize=GC_RE2_HEADER_LEN+sizeof(float);

  /* calculate read size for element data on each MPI rank */
  int readSize=nelt*elemDataSize;
  if(rank==0) readSize+=headerSize;

  char *buf=(char*)calloc(readSize,sizeof(char)),*buf0=buf;
  MPI_Status st;
  int err=MPI_File_read_ordered(file,buf,readSize,MPI_BYTE,&st);
  if(err) return 1;

  if(rank==0) buf0+=headerSize;

  /* initialize array */
  size_t nUnits=nelt*nVertex;
  exaArrayInit(&mesh->elements,struct Point_private,nUnits);
  exaArraySetSize(mesh->elements,nUnits);
  Point ptr=exaArrayGetPointer(mesh->elements);

  /* read elements for each rank */
  double x[GC_MAX_VERTICES],y[GC_MAX_VERTICES],\
    z[GC_MAX_VERTICES];
  int i,j,k;
  for(i=0;i<nelt;i++){
    // skip group id
    buf0+=sizeof(double);
    readT(x,buf0,double,nVertex); buf0+=sizeof(double)*nVertex;
    readT(y,buf0,double,nVertex); buf0+=sizeof(double)*nVertex;
    if(nDim==3){
      readT(z,buf0,double,nVertex);
      buf0+=sizeof(double)*nVertex;
    }

    for(k=0;k<nVertex;k++){
      j=PRE_TO_SYM_VERTEX[k];
      ptr->x[0]=x[j],ptr->x[1]=y[j];
      if(nDim==3)
        ptr->x[2]=z[j];
      ptr->elementId =start+i;
      ptr->sequenceId=nVertex*(start+i)+k;
      ptr->origin    =rank;
      ptr++;
    }
  }
  free(buf);

  return 0;
}

int readRe2Boundaries(exaHandle h,Mesh mesh,MPI_File file){
  exaInt rank=exaRank(h);
  exaInt size=exaSize(h);
  MPI_Comm comm=exaGetMPIComm(h);

  int nelt=mesh->nelt;
  int nelgt=mesh->nelgt;
  int nDim=mesh->nDim;
  int nVertex=mesh->nVertex;

  int elemDataSize=nVertex*nDim*sizeof(double)+sizeof(double);
  int headerSize=GC_RE2_HEADER_LEN+sizeof(float);

  MPI_Status st;
  char bufL[8];

  /* calculate offset for the curve side data */
  MPI_Offset curveOffset=headerSize+nelgt*elemDataSize;
  if(rank==0)
    MPI_File_read_at(file,curveOffset,bufL,sizeof(long),\
        MPI_BYTE,&st);
  MPI_Bcast(bufL,sizeof(long),MPI_BYTE,0,comm);

  double ncurvesD; readT(&ncurvesD,bufL,long,1);
  long ncurves=ncurvesD;

  /* calculate offset for boundary conditions data */
  MPI_Offset boundaryOffset=curveOffset+sizeof(long)+\
    sizeof(long)*8*ncurves;
  if(rank==0)
    MPI_File_read_at(file,boundaryOffset,bufL,sizeof(long),\
      MPI_BYTE,&st);
  MPI_Bcast(bufL,sizeof(long),MPI_BYTE,0,comm);

  double nbcsD; readT(&nbcsD,bufL,long,1); long nbcs=nbcsD;
  
  int nbcsLocal=nbcs/size,nrem=nbcs-nbcsLocal*size;
  nbcsLocal+=(rank<nrem ? 1 : 0);
  exaDebug(h,"rank=%d reads %d boundary faces.\n",rank,
      nbcsLocal);

  exaLong out[2][1],buff[2][1],in[1];
  in[0]=nbcsLocal;
  exaScan(h,out,in,buff,1,exaLong_t,exaAddOp);
  exaLong start=out[0][0];

  int offset=boundaryOffset+sizeof(long)+start*8*sizeof(long);
  int readSize=nbcsLocal*sizeof(long)*8;
  char *buf=calloc(readSize,sizeof(char)),*buf0=buf;
  MPI_File_read_at_all(file,offset,buf,readSize,MPI_BYTE,&st);

  exaArrayInit(&mesh->boundary,struct Boundary_private,0);

  double tmp[5];
  char cbc[4];
  struct Boundary_private boundary;
  exaInt i;
  for(i=0;i<nbcsLocal;i++){
    readT(tmp,buf0,long,1);buf0+=sizeof(long);
    boundary.elementId=tmp[0]-1;

    readT(tmp,buf0,long,1);buf0+=sizeof(long);
    boundary.faceId=PRE_TO_SYM_FACE[(long)tmp[0]-1];

    readT(tmp,buf0,long,5);buf0+=5*sizeof(long);
    readT(cbc,buf0,char,3);buf0+=sizeof(long);
    cbc[3]='\0';

    if(strcmp(cbc,GC_PERIODIC)==0){
      boundary.bc[0]=(long)tmp[0]-1;
      boundary.bc[1]=PRE_TO_SYM_FACE[(long)tmp[1]-1];
      exaArrayAppend(mesh->boundary,&boundary);
    }
  }
  free(buf);
}

int genConReadRe2File(exaHandle h,Mesh *mesh_,char *fileName){
  int nelt,nDim,nVertex;
  int errs;

  exaInt rank=exaRank(h);
  exaInt size=exaSize(h);
  MPI_Comm comm=exaGetMPIComm(h);

  MPI_File file;
  int err=MPI_File_open(comm,fileName,MPI_MODE_RDONLY,\
    MPI_INFO_NULL,&file);
  if(err){
    if(rank==0)
      printf("%s:%d Error opening file: %s\n",
        __FILE__,__LINE__,fileName);
    errs++;
    MPI_Abort(comm,911);
  }

  exaMalloc(1,mesh_);
  Mesh mesh=*mesh_;

  readRe2Header(h,mesh,file);
  readRe2Coordinates(h,mesh,file);
  readRe2Boundaries(h,mesh,file);
  transferBoundaryFaces(h,mesh);

  err=MPI_File_close(&file);
  if(err) errs++;

  MPI_Barrier(comm);

  return errs;
}

int genConWriteCo2File(exaHandle h,Mesh mesh,char *fileName){
  const char version[5]="#v001";
  const float test=6.54321;

  exaInt rank=exaRank(h);
  exaInt size=exaSize(h);
  MPI_Comm comm=exaGetMPIComm(h);

  int nVertex=mesh->nVertex;
  int nDim=mesh->nDim;
  exaInt nelt=mesh->nelt;
  exaLong nelgt=mesh->nelgt;
  exaLong nelgv=mesh->nelgv;

  int errs=0;

  MPI_File file;
  int err=MPI_File_open(comm,fileName,\
    MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&file);
  if(err){
    if(rank==0)
      printf("%s:%d Error opening file: %s for writing.\n",
        __FILE__,__LINE__,fileName);
    errs++;
    MPI_Abort(comm,911);
  }

  exaLong out[2][1],buff[2][1],in[1];
  in[0]=nelt;
  exaScan(h,out,in,buff,1,exaLong_t,exaAddOp);
  exaLong start=out[0][0];

  int writeSize=nelt*(nVertex+1)*sizeof(int);
  int headerSize=GC_CO2_HEADER_LEN+sizeof(float);
  if(rank==0) writeSize+=headerSize;

  char *buf=(char*)calloc(writeSize,sizeof(char)),*buf0=buf;
  MPI_Status st;
  if(rank==0){
    sprintf(buf0,"%5s%12d%12d%12d",version,(int)nelgt,\
        (int)nelgv,nVertex);
    exaDebug(h,"%5s%12d%12d%12d\n",version,(int)nelgt,\
        (int)nelgv,nVertex);
    memset(buf0+strlen(buf0),' ',GC_CO2_HEADER_LEN-strlen(buf0));
    buf0[GC_CO2_HEADER_LEN]='\0';
    buf0+=GC_CO2_HEADER_LEN;
    memcpy(buf0,&test,sizeof(float)),buf0+=sizeof(float);
  }

  Point ptr=exaArrayGetPointer(mesh->elements);
  int i,k,temp;
  for(i=0;i<nelt;i++){
    temp=ptr->elementId+1;
    writeInt(buf0,temp); buf0+=sizeof(int);
    exaDebug(h,"%lld",temp);
    for(k=0;k<nVertex;k++){
      temp=ptr->globalId+1;
      writeInt(buf0,temp); buf0+=sizeof(int);
      exaDebug(h," %lld",temp);
      ptr++;
    }
    exaDebug(h,"\n");
  }

  err=MPI_File_write_ordered(file,buf,writeSize,MPI_BYTE,&st);
  if(err) errs++;

  err=MPI_File_close(&file);
  if(err) errs++;

  MPI_Barrier(comm);
  free(buf);

  return errs;
}

#undef readT
#undef writeInt
