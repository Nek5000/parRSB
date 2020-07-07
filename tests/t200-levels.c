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
  mgData d; mgSetup(c,M,&d); uint nlevels=d->nLevels;

  GenmapScalar *x=d->x,*y=d->y,*buf=d->buf;

  for(i=0; i<d->level_off[nlevels]; i++)
    y[i]=10.,x[i]=1.0;

  for(i=0; i<nlevels; i++)
    parMatApply(y+d->level_off[i],d->levels[i]->M,x+d->level_off[i],buf);

  for(j=0; j<d->level_off[nlevels]; j++)
    assert(fabs(y[j])<GENMAP_TOL);

  mgFree(d);

  GenmapFinalize(gh);
  MeshFree(mesh);

  comm_free(&comm);
  MPI_Finalize();

  return 0;
}
