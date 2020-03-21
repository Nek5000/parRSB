/*
Parition mesh using Nek5000's geometry (re2) file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "exa-impl.h"
#include "gencon-impl.h"

#include "quality.h"
#include "parRSB.h"

#define CHECK_ERR(ierr) do{\
  if(ierr){\
    MPI_Finalize();\
    return ierr;\
  }\
} while(0);

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int id; MPI_Comm_rank(MPI_COMM_WORLD,&id);
  int color=0;

  if (argc!=3){
    if(id==0) printf("usage: ./example <#nread> <re2 file>\n");
    return EXIT_FAILURE;
  } 

  int nRead=atoi(argv[1]);
  char* reaFile=argv[2];

  if(id<nRead) color=1;
  MPI_Comm commRead;
  MPI_Comm_split(MPI_COMM_WORLD,color,id,&commRead);

  exaHandle h;
  exaInit(&h,commRead,"/host");

  int ierr=0;
  Mesh mesh;
  if(color==1)
    ierr=genConReadRe2File(h,&mesh,reaFile);
  CHECK_ERR(ierr);

  int nelt=mesh->nelt;
  int ndim=mesh->nDim;
  int nv  =mesh->nVertex;

  int *part;   exaMalloc(nelt,&part);
  double *vtx; exaCalloc(nelt*ndim,&vtx);

  Point points=exaArrayGetPointer(mesh->elements);
  int e,n,d;
  for(e=0;e<nelt;e++)
    for(n=0;n<nv;n++)
      for(d=0;d<ndim;d++)
        vtx[e*ndim+d]+=points[e*nv+n].x[d];

  int options[3];
  options[0] = 1; /* use custom options */
  options[1] = 3; /* debug level        */
  options[2] = 0; /* not used           */

  ierr=parRCB_partMesh(part,vtx,nelt,ndim,options,commRead);
  CHECK_ERR(ierr);

  /* redistribute data */
  for(e=0;e<nelt;e++)
    for(n=0;n<nv;n++)
      points[e*nv+n].proc=part[e];

  exaArrayTransfer(mesh->elements,
    offsetof(struct Point_private,proc),0,exaGetComm(h));
  exaSortArray(mesh->elements,exaLong_t,
    offsetof(struct Point_private,sequenceId));

  exaFree(vtx);
  exaFree(part);

  return 0;
}
