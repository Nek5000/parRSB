#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <exa.h>

#include <genmap-impl.h>
#include <gencon-impl.h>

int main(int argc,char *argv[]){
  MPI_Init(&argc,&argv);

  exaHandle h; exaInit(&h,MPI_COMM_WORLD,"/host");
  int rank=exaRank(h);
  int size=exaSize(h);

  if(argc!=2){
    if(rank==0) printf("Usage: ./%s <co2 file>\n",argv[0]);
    MPI_Finalize();
    exit(1);
  }

  Mesh mesh;
  readCo2File(h,&mesh,argv[1]);

  GenmapHandle gh; GenmapInit(&gh,MPI_COMM_WORLD);
  GenmapSetNLocalElements(gh,mesh->nelt);
  GenmapScan(gh,GenmapGetGlobalComm(gh));
  GenmapSetNVertices(gh,mesh->nVertex);

  GenmapElements e=GenmapGetElements(gh);
  GenmapLong start=GenmapGetLocalStartIndex(gh);
  GenmapInt id=GenmapCommRank(GenmapGetGlobalComm(gh));

#if 0
  Element elements=MeshGetElements(m);
  GenmapInt i,j;
  for(i=0;i<nelt;i++){
    for(j=0;j<nv;j++)
      e[i].vertices[j]=elements[count].globalId;
  }
  MeshFree(m);
#endif

  GenmapVector weights; GenmapCreateVector(&weights,mesh->nelt);

  GenmapInitLaplacian(gh,GenmapGetGlobalComm(gh),weights);

  GenmapDestroyVector(weights);
  GenmapFinalize(gh);

  exaFinalize(h);
  MPI_Finalize();

  return 0;
}
