#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <exa.h>

#include <genmap-impl.h>
#include <gencon-impl.h>
#include <genmap-multigrid-precon.h>

int main(int argc,char *argv[]){
  MPI_Init(&argc,&argv);
  struct comm comm; comm_init(&comm,MPI_COMM_WORLD);
  int rank=comm.id,size=comm.np;

  if(argc!=2){
    if(rank==0) printf("Usage: ./%s <co2 file>\n",argv[0]);
    MPI_Finalize();
    exit(1);
  }

  Mesh mesh;
  readCo2File(&mesh,argv[1],&comm);

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
  GenmapInitLaplacian(gh,c);
  GenmapLaplacian(gh,c,u,v);

  for(i=0;i<mesh->nelt;i++)
    assert(fabs(v->data[i])<GENMAP_TOL);

  GenmapDestroyVector(weights);

  for(i=0;i<mesh->nelt;i++)
    v->data[i]=100; // Random number to reset v for next test

  /* Test Laplacian based on CSR representation */
  csr_mat M; csr_mat_setup(gh,c,&M);

  GenmapScalar *x,*y,*buf;
  GenmapMalloc(mesh->nelt,&x); GenmapMalloc(mesh->nelt,&y);
  GenmapUInt nnz=M->row_off[M->rn]; GenmapMalloc(nnz,&buf);

  for(i=0; i<mesh->nelt; i++) x[i]=1.0;

  buffer bfr; buffer_init(&bfr,1024);
  csr_mat_gather(M,M->gsh,x,buf,&bfr);
  csr_mat_apply(y,M,buf);
  buffer_free(&bfr);

  for(i=0;i<mesh->nelt;i++)
    assert(fabs(y[i])<GENMAP_TOL);

  GenmapFree(buf); GenmapFree(x); GenmapFree(y);
  csr_mat_free(M);

  GenmapDestroyVector(v); GenmapDestroyVector(u);

  GenmapFinalize(gh);
  MeshFree(mesh);

  comm_free(&comm);
  MPI_Finalize();

  return 0;
}
