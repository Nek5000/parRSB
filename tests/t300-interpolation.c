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

  /* Setup mesh */
  GenmapElements e=GenmapGetElements(gh);
  Element me      =MeshGetElements(mesh);
  GenmapInt i,j;
  for(i=0;i<mesh->nelt;i++)
    for(j=0;j<mesh->nVertex;j++)
      e[i].vertices[j]=me[i].vertex[j].globalId;

  /* Setup CSR on fine level */
  GenmapComm c=GenmapGetGlobalComm(gh);
  parMat M; parMatSetup(gh,c,&M);

  /* Setup MG levels */
  mgData d; mgSetup(c,M,&d);

  uint nlevels=d->nLevels; mgLevel *l=d->levels;
  uint *lvl_off=d->level_off; GenmapScalar *x=d->x;

  for(i=lvl_off[0]; i<lvl_off[1]; i++)
    x[i]=1.0;

  buffer buf; buffer_init(&buf,1024);
  GenmapScalar v; GenmapULong rn=M->rn;
  for(i=0; i<nlevels-1; i++){
    if(lvl_off[i]==lvl_off[i+1])
      continue;

    gs(x+lvl_off[i],gs_double,gs_add,1,l[i]->J,&buf);

    v=0.0;
    for(j=lvl_off[i+1]; j<lvl_off[i+2]; j++)
      v+=x[j];

    printf("x: ");
    for(j=lvl_off[i]; j<lvl_off[i+2]; j++)
      printf(" %lf",x[j]);
    printf("\nv=%lf rn=%d\n",v,rn);

    assert(fabs(v-rn)<GENMAP_TOL);

    for(j=lvl_off[i+1]; j<lvl_off[i+2]; j++)
      x[j]=1.0;
  }

  buffer_free(&buf);

  mgFree(d);

  GenmapFinalize(gh);

  MeshFree(mesh);
  exaFinalize(h);

  MPI_Finalize();

  return 0;
}
