#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "genmap-impl.h"
#include "conReader.h"

int main(int argc,char *argv[]){
  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  int size; MPI_Comm_size(MPI_COMM_WORLD,&size);

  if(argc!=2){
    if(rank==0) printf("Usage: ./%s <co2 file>\n",argv[0]);
    MPI_Finalize();
    exit(1);
  }

  char *conFile=argv[1];
  struct con con;
  conRead(conFile,&con,MPI_COMM_WORLD);

  GenmapInt nel=con.nel;
  int nv=con.nv;

  GenmapHandle h; GenmapInit(&h,MPI_COMM_WORLD);
  GenmapSetNLocalElements(h,nel);
  GenmapScan(h,GenmapGetGlobalComm(h));
  GenmapSetNVertices(h,nv);

  GenmapElements e=GenmapGetElements(h);
  GenmapLong start=GenmapGetLocalStartIndex(h);
  GenmapInt id=GenmapCommRank(GenmapGetGlobalComm(h));

  GenmapInt i,j;
  for(i=0;i<nel;i++){
    e[i].origin=id;
    for(j=0;j<nv;j++)
      e[i].vertices[j]=con.vl[i*nv+j];
  }

  GenmapVector weights; GenmapCreateVector(&weights,nel);

  GenmapInitLaplacian(h,GenmapGetGlobalComm(h),weights);

  GenmapDestroyVector(weights);
  GenmapFinalize(h);

  conFree(&con);

  MPI_Finalize();

  return 0;
}
