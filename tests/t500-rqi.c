#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <genmap-impl.h>
#include <gencon-impl.h>
#include <genmap-multigrid-precon.h>
#include <parRSB.h>

int main(int argc,char *argv[]){
  MPI_Init(&argc,&argv);
  struct comm comm; comm_init(&comm,MPI_COMM_WORLD);
  int rank=comm.id,size=comm.np;

  exaHandle h; exaInit(&h,MPI_COMM_WORLD,"/host");

  if(argc>6){
    if(rank==0)
      printf("Usage: ./%s <re2 file> <co2 file> 1 2 3 \n",argv[0]);
    MPI_Finalize();
    exit(1);
  }

  Mesh mesh;
  readRe2File(h,&mesh,argv[1]);
  read_co2_file(mesh,argv[2],&comm);

  GenmapInt i,j;

  int seq_g=(argc>3)?atoi(argv[3]):1;
  int rcb_g=(argc>4)?atoi(argv[4]):1;
  int rcb_l=(argc>5)?atoi(argv[5]):0;

  if(seq_g){
    exaSort(mesh->elements,
      exaULong_t,offsetof(struct Point_private,sequenceId),
      exaSortAlgoBinSort,1,exaGetComm(h));
    mesh->nelt=exaArrayGetSize(mesh->elements)/mesh->nVertex;
  }

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

  me=(Point)MeshGetElements(mesh);

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
#if 0
      printf("centroid: %lf %lf %lf\n",\
          ptr[i].coord[0],ptr[i].coord[1],ptr[i].coord[2]);
#endif
    }

    buffer buf; buffer_init(&buf,1024);
    uint s1=0,e1=mesh->nelt;
    rcb_local(&a,s1,e1,mesh->nDim,&buf);
    ptr=a.ptr;

    for(i=0; i<mesh->nelt; i++){
      ptr[i].proc=i;
#if 0
      printf("i=%d %lf %lf %lf\n",i,ptr[i].coord[0],\
          ptr[i].coord[1],ptr[i].coord[2]);
#endif
    }

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

  me=(Point)MeshGetElements(mesh);

  GenmapHandle gh; GenmapInit(&gh,MPI_COMM_WORLD);
  GenmapSetNLocalElements(gh,mesh->nelt);
  GenmapSetNVertices(gh,mesh->nVertex);

  /* Setup mesh */
  GenmapElements e=GenmapGetElements(gh);
  me=(Point)MeshGetElements(mesh);
  for(i=0;i<mesh->nelt;i++)
    for(j=0;j<mesh->nVertex;j++)
      e[i].vertices[j]=me[i*mesh->nVertex+j].globalId;

  GenmapComm c=GenmapGetGlobalComm(gh);
  GenmapInitLaplacian(gh,c);

  GenmapVector x; GenmapCreateVector(&x,mesh->nelt);
  GenmapVector r; GenmapCreateVector(&r,mesh->nelt);

  srand(time(0));
  for(i=0; i<mesh->nelt; i++){
#if 0
    x->data[i]=me[i*mesh->nVertex].elementId;
    printf(" %lf",x->data[i]);
#else
    x->data[i]=rand()%100/50.;
#endif
  }

  //csr_mat_print(c->M,&c->gsc);

  GenmapLong nelg=GenmapGetNGlobalElements(gh);
  GenmapOrthogonalizebyOneVector(gh,c,x,nelg);

  GenmapScalar norm=GenmapDotVector(x,x);
  GenmapGop(c,&norm,1,GENMAP_SCALAR,GENMAP_SUM);
  //printf("norm(z)=%lf\n",sqrt(norm));

  GenmapScalar normi=1.0/sqrt(norm);
  GenmapAxpbyVector(x,x,0.0,x,normi);

  norm=GenmapDotVector(x,x);
  GenmapGop(c,&norm,1,GENMAP_SCALAR,GENMAP_SUM);
  //printf("norm(z)=%lf\n",sqrt(norm));

  rqi(gh,c,x,30,1,r);

  GenmapDestroyVector(x); GenmapDestroyVector(r);

  GenmapFinalize(gh);
  MeshFree(mesh);

  comm_free(&comm);
  exaFree(h); MPI_Finalize();

  return 0;
}