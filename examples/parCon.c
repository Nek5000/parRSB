#include <stdio.h>
#include <stdlib.h>

#include <gencon-impl.h>

int main(int argc,char *argv[]){
  MPI_Init(&argc,&argv);

  exaHandle h;
  exaInit(&h,MPI_COMM_WORLD,"/host");
  struct comm comm; comm_init(&comm,MPI_COMM_WORLD);

  if(argc<2) {
    if(exaRank(h)==0) printf("Usage: ./%s foo.re2 [tol]\n",argv[0]);
    exit(1);
  }

  Mesh mesh;
  read_re2_mesh(&mesh,argv[1],&comm);

  findMinNeighborDistance(h,mesh);

  double tol=0.2;
  if(argc==3) tol=atof(argv[2]);
  findSegments(h,mesh,tol);

  setGlobalID(h,mesh);
  sendBack(h,mesh);
  matchPeriodicFaces(h,mesh);

  char co2FileName[BUFSIZ]; strcpy(co2FileName,argv[1]);
  int len=strlen(co2FileName); assert(len>4);
  co2FileName[len-2]='o',co2FileName[len-3]='c';

  write_co2_file(mesh,co2FileName,&comm);

  MeshFree(mesh);

  comm_free(&comm);
  exaFinalize(h);

  MPI_Finalize();

  return 0;
}
