/*
Parition mesh using Nek5000's geometry (re2) file.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "exa-impl.h"
#include "gencon-impl.h"

#include "parRSB.h"

#define CHECK_ERR(ierr) do{\
  if(ierr){\
    MPI_Finalize();\
    return ierr;\
  }\
} while(0);

int main(int argc,char *argv[]){
  MPI_Init(&argc,&argv);

  int id; MPI_Comm_rank(MPI_COMM_WORLD,&id);
  int color=0;

  if(argc!=3){
    if(id==0) printf("usage: ./example <#nread> <re2 file>\n");
    return EXIT_FAILURE;
  } 

  int nRead=atoi(argv[1]);
  char* reaFile=argv[2];

  if(id<nRead) color=1;
  MPI_Comm commRead;
  MPI_Comm_split(MPI_COMM_WORLD,color,id,&commRead);

  if(color==1){
    exaHandle h;
    exaInit(&h,commRead,"/host");
    comm comm; comm_init(&comm,commRead);

    Mesh mesh;
    int ierr=read_re2_mesh(&mesh,reaFile,&comm);
    CHECK_ERR(ierr);

    int nelt=mesh->nelt;
    int ndim=mesh->nDim;
    int nv  =mesh->nVertex;

    int *part;   exaMalloc(nelt,&part);
    double *vtx; exaCalloc(nelt*ndim,&vtx);

    Point points=exaArrayGetPointer(mesh->elements);
    assert(exaArrayGetSize(mesh->elements)==nelt*nv);

    int e,n,d;
    for(e=0;e<nelt;e++){
      for(n=0;n<nv;n++){
        for(d=0;d<ndim;d++)
          vtx[e*ndim+d]+=points[e*nv+n].x[d];
      }
      for(d=0;d<ndim;d++) vtx[e*ndim+d]/=nv;
    }

    int options[3];
    options[0] = 1; /* use custom options */
    options[1] = 3; /* debug level        */
    options[2] = 0; /* not used           */

    ierr=parRCB_partMesh(part,vtx,nelt,nv,options,commRead);
    CHECK_ERR(ierr);

    /* Redistribute data */
    for(e=0;e<nelt;e++)
      for(n=0;n<nv;n++)
        points[e*nv+n].proc=part[e];

    exaArrayTransfer(mesh->elements,
      offsetof(struct Point_private,proc),0,exaGetComm(h));
    exaSortArray(mesh->elements,exaULong_t,
      offsetof(struct Point_private,sequenceId));

    /* Write output */
    int nPoints=exaArrayGetSize(mesh->elements);
    points=exaArrayGetPointer(mesh->elements);

    MPI_File file;
    int err=MPI_File_open(commRead,"out.part",\
      MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&file);

    int writeSize=nPoints*(ndim*sizeof(double)+sizeof(int));

    char header[BUFSIZ];
    sprintf(header,"%lld %d %d %ld %ld",mesh->nelgt,ndim,nv,
      sizeof(double),sizeof(int));

    int rank; MPI_Comm_rank(commRead,&rank);
    if(rank==0) writeSize+=128;

    char *buf,*bufPtr; exaCalloc(writeSize,&buf),bufPtr=buf;

    if(rank==0){
      strcpy(bufPtr,header);
      memset(bufPtr+strlen(header),' ',128-strlen(header));
      bufPtr+=128;
    }

    for(e=0;e<nPoints;e++){
      memcpy(bufPtr,points[e].x,sizeof(double)*ndim);
      bufPtr+=ndim*sizeof(double);
      memcpy(bufPtr,&rank,sizeof(int));
      bufPtr+=sizeof(int);
    }

    MPI_Status st;
    err|=MPI_File_write_ordered(file,buf,writeSize,MPI_BYTE,&st);
    err|=MPI_File_close(&file);

    MPI_Barrier(commRead);

    assert(err==0);

    MeshFree(mesh);

    exaFree(buf);
    exaFree(vtx);
    exaFree(part);
    comm_free(&comm);
    exaFinalize(h);
  }

  MPI_Finalize();

  return 0;
}
