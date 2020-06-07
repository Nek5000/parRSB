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
  GenmapSetNVertices(gh,mesh->nVertex);

  GenmapElements e=GenmapGetElements(gh);
  Element me      =MeshGetElements(mesh);
  GenmapInt i,j;
  for(i=0;i<mesh->nelt;i++){
    for(j=0;j<mesh->nVertex;j++)
      e[i].vertices[j]=me[i].vertex[j].globalId;
  }

  GenmapVector weights,u,v;
  GenmapCreateVector(&weights,mesh->nelt);
  GenmapCreateVector(&u      ,mesh->nelt);
  GenmapCreateVector(&v      ,mesh->nelt);

  for(i=0;i<mesh->nelt;i++)
    u->data[i]=1.0;

  GenmapInitLaplacian(gh,GenmapGetGlobalComm(gh),weights);
  GenmapLaplacian(gh,GenmapGetGlobalComm(gh),u,weights,v);

  printf("v: ");
  for(i=0;i<mesh->nelt;i++)
    printf(" %lf",v->data[i]);
  printf("\n");

  GenmapDestroyVector(weights);
  GenmapDestroyVector(v      );
  GenmapDestroyVector(u      );

  GenmapFinalize(gh);
  MeshFree(mesh);
  exaFinalize(h);

  MPI_Finalize();

  return 0;
}
