#include <exasort.h>
#include <exa-memory.h>
#include <time.h>

typedef struct{
  exaScalar ds;
  uint proc;
} Data;

#define N 10

int main(int argc,char *argv[]){
  MPI_Init(&argc,&argv);

  exaHandle h;
  exaInit(&h,MPI_COMM_WORLD,"/host");

  srand(time(0));

  exaArray arr; exaArrayInit(&arr,Data,N);

  int i; Data d;
  for(i=0; i<N; i++)
      d.ds=(rand()%100)/100.0,exaArrayAppend(arr,&d);

  exaSort(arr,exaScalar_t,offsetof(Data,ds),exaSortAlgoBinSort,1,
      exaGetComm(h));

  assert(exaArrayGetSize(arr)==N);

  exaDestroy(arr);
  exaFinalize(h);
  MPI_Finalize();

  return 0;
}
