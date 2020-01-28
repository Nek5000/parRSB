#include <stdio.h>
#include <stdlib.h>

#include <gencon-impl.h>

int main(int argc,char *argv[]){
  MPI_Init(&argc,&argv);

  exaHandle h;
  exaInit(&h,MPI_COMM_WORLD,"/host");

  if(argc!=2) {
    if(exaRank(h)==0) printf("Usage: ./example foo.re2\n");
    exit(0);
  }

  Mesh mesh;
  genConReadRe2File(h,&mesh,argv[1]);

  findMinNeighborDistance(h,mesh);
  findSegments(h,mesh,1.e-2);
  setGlobalID(h,mesh);
  sendBack(h,mesh);
  matchPeriodicFaces(h,mesh);

  char co2FileName[BUFSIZ]; strcpy(co2FileName,argv[1]);
  int len=strlen(co2FileName); assert(len>4);
  co2FileName[len-2]='o',co2FileName[len-3]='c';

  genConWriteCo2File(h,mesh,co2FileName);

  freeMesh(mesh);

  exaFinalize(h);

  MPI_Finalize();

  return 0;
}
