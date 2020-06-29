#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <exa.h>

#include <genmap-impl.h>
#include <gencon-impl.h>
#include <genmap-multigrid-precon.h>

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

  /* Test GenmapLaplacian based on gather-scatter */
  GenmapElements e=GenmapGetElements(gh);
  Element me      =MeshGetElements(mesh);
  GenmapInt i,j;
  for(i=0;i<mesh->nelt;i++)
    for(j=0;j<mesh->nVertex;j++)
      e[i].vertices[j]=me[i].vertex[j].globalId;

  GenmapVector weights,u,v;
  GenmapCreateVector(&weights,mesh->nelt);
  GenmapCreateVector(&u      ,mesh->nelt);
  GenmapCreateVector(&v      ,mesh->nelt);

  for(i=0;i<mesh->nelt;i++)
    u->data[i]=1.0;

  GenmapComm c=GenmapGetGlobalComm(gh);
  GenmapInitLaplacian(gh,c,weights);
  GenmapLaplacian(gh,c,u,weights,v);

  for(i=0;i<mesh->nelt;i++)
    assert(fabs(v->data[i])<GENMAP_TOL);

  GenmapDestroyVector(weights);

  for(i=0;i<mesh->nelt;i++)
    v->data[i]=100; // Random number to reset v for next test

#if 1
  /* Test Laplacian based on CSR representation */
  parMat M; parMatSetup(gh,c,&M);

  GenmapScalar *x,*y,*buf;
  GenmapMalloc(mesh->nelt,&x); GenmapMalloc(mesh->nelt,&y);
  GenmapUInt nnz=M->rowOffsets[M->rn]; GenmapMalloc(nnz,&buf);

  for(i=0; i<mesh->nelt; i++) x[i]=1.0;

  parMatApply(y,M,x,buf);
  for(i=0;i<mesh->nelt;i++)
    assert(fabs(y[i])<GENMAP_TOL);

  mgData d;
  mgSetup(c,M,&d);

  parMatFree(M);

  GenmapFree(buf); GenmapFree(x); GenmapFree(y);
#endif

  GenmapDestroyVector(v); GenmapDestroyVector(u);

  GenmapFinalize(gh);
  MeshFree(mesh);
  exaFinalize(h);

  MPI_Finalize();

  return 0;
}
