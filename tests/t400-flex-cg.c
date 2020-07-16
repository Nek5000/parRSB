#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <genmap-impl.h>
#include <gencon-impl.h>
#include <genmap-multigrid-precon.h>

int main(int argc,char *argv[]){
  MPI_Init(&argc,&argv);
  struct comm comm; comm_init(&comm,MPI_COMM_WORLD);
  int rank=comm.id,size=comm.np;

  exaHandle h; exaInit(&h,MPI_COMM_WORLD,"/host");

  if(argc!=3){
    if(rank==0) printf("Usage: ./%s <re2 file> <co2 file>\n",argv[0]);
    MPI_Finalize();
    exit(1);
  }

  Mesh mesh;
  readRe2File(h,&mesh,argv[1]);
  read_co2_file(mesh,argv[2],&comm);

  GenmapHandle gh; GenmapInit(&gh,MPI_COMM_WORLD);

  GenmapSetNLocalElements(gh,mesh->nelt);
  GenmapSetNVertices(gh,mesh->nVertex);

  /* Setup mesh */
  GenmapElements e=GenmapGetElements(gh);
  Point   me      =(Point) MeshGetElements(mesh);
  GenmapInt i,j;
  for(i=0;i<mesh->nelt;i++)
    for(j=0;j<mesh->nVertex;j++)
      e[i].vertices[j]=me[i*mesh->nVertex+j].globalId;

  /* Setup CSR on fine level */
  GenmapComm c=GenmapGetGlobalComm(gh);
  csr_mat M; csr_mat_setup(gh,c,&M);

  /* Setup MG levels */
  mgData d; mgSetup(c,M,&d);

  GenmapVector r,x,x0,weights;
  GenmapCreateVector(&weights,mesh->nelt);
  GenmapCreateVector(&r      ,mesh->nelt);
  GenmapCreateVector(&x      ,mesh->nelt);
  GenmapCreateVector(&x0     ,mesh->nelt);

  srand(time(0));
  for(i=0; i<mesh->nelt; i++)
    x->data[i]=rand()%100/100.0,x0->data[i]=0.0;

  GenmapLong nelg=GenmapGetNGlobalElements(gh);
  GenmapOrthogonalizebyOneVector(gh,c,x,nelg);

  GenmapInitLaplacian(gh,c,weights);
  GenmapLaplacian(gh,c,x,weights,r);

  i=flex_cg(gh,c,d,r,weights,100,x0);
  if(rank==0)
    printf("Flex-CG iterations: %d\n",i);

  for(i=0; i<mesh->nelt; i++){
    GenmapScalar e=x->data[i]-x0->data[i];
    assert(fabs(e)<1e-10);
  }

  GenmapDestroyVector(r ); GenmapDestroyVector(x);
  GenmapDestroyVector(x0); GenmapDestroyVector(weights);

  mgFree(d);

  GenmapFinalize(gh);
  MeshFree(mesh);

  comm_free(&comm);
  exaFree(h);
  MPI_Finalize();

  return 0;
}
