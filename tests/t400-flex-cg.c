#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <genmap-impl.h>
#include <gencon-impl.h>
#include <genmap-multigrid-precon.h>
#include <parRSB.h>

int main(int argc,char *argv[]){
  MPI_Init(&argc,&argv);

  struct comm comm;
  comm_init(&comm,MPI_COMM_WORLD);
  int rank=comm.id,size=comm.np;

  exaHandle h; exaInit(&h,MPI_COMM_WORLD,"/host");

  if(argc>5){
    if(rank==0)
      printf("Usage: ./%s <re2 file> <co2 file> global local\n",argv[0]);
    MPI_Finalize();
    exit(1);
  }

  Mesh mesh;
  read_re2_mesh(&mesh,argv[1],&comm);
  read_co2_file( mesh,argv[2],&comm);

  GenmapInt i,j;

  int rcb_g=(argc>3)?atoi(argv[3]):1;
  int rcb_l=(argc>4)?atoi(argv[4]):1;

  //partition
  int      *part; GenmapMalloc(mesh->nelt              ,&part  );
  uint    *upart; GenmapMalloc(mesh->nelt*mesh->nVertex,&upart );
  double *coords; GenmapMalloc(mesh->nelt*mesh->nDim   ,&coords);

  int nDim=mesh->nDim;
  Point me=(Point)MeshGetElements(mesh);

  if(rcb_g){
    for(i=0; i<mesh->nelt; i++){
      coords[i*nDim+0]=0.0;
      coords[i*nDim+1]=0.0;
      if(nDim==3)
        coords[i*nDim+2]=0.0;
      for(j=0; j<mesh->nVertex; j++){
        coords[i*nDim+0]+=me[i*mesh->nVertex+j].x[0];
        coords[i*nDim+1]+=me[i*mesh->nVertex+j].x[1];
        if(nDim==3)
          coords[i*nDim+2]+=me[i*mesh->nVertex+j].x[2];
      }
      coords[i*nDim+0]/=mesh->nVertex;
      coords[i*nDim+1]/=mesh->nVertex;
      if(nDim==3)
        coords[i*nDim+2]/=mesh->nVertex;
    }

    int options[3]; options[0]=options[1]=options[2]=0;
    parRCB_partMesh(part,coords,mesh->nelt,mesh->nVertex,options,
      MPI_COMM_WORLD);

    for(i=0; i<mesh->nelt; i++)
      for(j=0; j<mesh->nVertex; j++)
        upart[i*mesh->nVertex+j]=part[i];

    exaArrayTransferExt(mesh->elements,upart,exaGetComm(h));
    exaSortArray(mesh->elements,exaULong_t,
      offsetof(struct Point_private,sequenceId));
    mesh->nelt=exaArrayGetSize(mesh->elements)/mesh->nVertex;
  }

  if(rcb_l){
    struct array a; array_init(elm_rcb,&a,mesh->nelt); a.n=mesh->nelt;
    elm_rcb *ptr=a.ptr;
    for(i=0; i<mesh->nelt; i++){
      ptr[i].orig=i;
      ptr[i].coord[0]=0.0;
      ptr[i].coord[1]=0.0;
      if(nDim==3)
        ptr[i].coord[2]=0.0;
      for(j=0; j<mesh->nVertex; j++){
        ptr[i].coord[0]+=me[i*mesh->nVertex+j].x[0];
        ptr[i].coord[1]+=me[i*mesh->nVertex+j].x[1];
        if(nDim==3)
          ptr[i].coord[2]+=me[i*mesh->nVertex+j].x[2];
      }
      ptr[i].coord[0]/=mesh->nVertex;
      ptr[i].coord[1]/=mesh->nVertex;
      if(nDim==3)
        ptr[i].coord[2]/=mesh->nVertex;
    }

    buffer buf; buffer_init(&buf,1024);
    uint s1=0,e1=mesh->nelt;
    rcb_local(&a,s1,e1,mesh->nDim,&buf);
    ptr=a.ptr;

    for(i=0; i<mesh->nelt; i++) ptr[i].proc=i;

    sarray_sort(elm_rcb,a.ptr,a.n,orig,0,&buf);
    ptr=a.ptr;

    Point pp=exaArrayGetPointer(mesh->elements);
    int cnt=0;
    for(i=0; i<mesh->nelt; i++)
      for(j=0; j<mesh->nVertex; j++){
        pp[cnt].proc=ptr[i].proc;
        cnt++;
      }

    buffer_free(&buf);
    array_free(&a);

    exaSortArray2(mesh->elements,
        exaUInt_t ,offsetof(struct Point_private,proc      ),
        exaULong_t,offsetof(struct Point_private,sequenceId));
  }

  free(upart);
  free(part);
  free(coords);

  GenmapHandle gh; GenmapInit(&gh,MPI_COMM_WORLD);
  GenmapSetNLocalElements(gh,mesh->nelt);
  GenmapSetNVertices(gh,mesh->nVertex);

  /* Setup mesh */
  GenmapElements e=GenmapGetElements(gh);
  me=(Point)MeshGetElements(mesh);
  for(i=0;i<mesh->nelt;i++)
    for(j=0;j<mesh->nVertex;j++)
      e[i].vertices[j]=me[i*mesh->nVertex+j].globalId;

  /* Setup CSR on fine level */
  GenmapComm c=GenmapGetGlobalComm(gh);
  struct array *entries=GenmapFindNeighbors(gh,c);
  csr_mat M; csr_mat_setup(entries,&c->gsc,&M);
  array_free(entries); free(entries);

  /* Setup MG levels */
  mgData d; mgSetup(c,M,&d);

  GenmapVector r,x,x0;
  GenmapCreateVector(&r      ,mesh->nelt);
  GenmapCreateVector(&x      ,mesh->nelt);
  GenmapCreateVector(&x0     ,mesh->nelt);

  srand(time(0));
  for(i=0; i<mesh->nelt; i++)
    x->data[i]=rand()%100/50.,x0->data[i]=0.0;

  GenmapLong nelg=GenmapGetNGlobalElements(gh);
  GenmapOrthogonalizebyOneVector(gh,c,x,nelg);

  GenmapInitLaplacian(gh,c);
  GenmapLaplacian(gh,c,x,r);

  i=flex_cg(gh,c,d,r,30,1,x0);
  if(rank==0)
    printf("Flex-CG iterations: %d\n",i);

  for(i=0; i<mesh->nelt; i++){
    GenmapScalar e=x->data[i]-x0->data[i];
    assert(fabs(e)<1e-10);
  }

  GenmapDestroyVector(r );
  GenmapDestroyVector(x );
  GenmapDestroyVector(x0);

  mgFree(d);

  GenmapFinalize(gh);
  MeshFree(mesh);

  comm_free(&comm);
  exaFree(h);
  MPI_Finalize();

  return 0;
}
