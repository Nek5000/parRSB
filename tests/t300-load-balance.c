#include <genmap-impl.h>
#include <sort-impl.h>
#include <time.h>

typedef struct{
  double ds;
  uint proc;
} Data;

#define N 10

int main(int argc,char *argv[]){
  MPI_Init(&argc,&argv);
  struct comm c; comm_init(&c,MPI_COMM_WORLD);

  srand(time(0));

  struct array arr; array_init(Data,&arr,N);

  int i,cnt; Data d; Data *ptr=arr.ptr;
  for(i=cnt=0; i<N; i++)
      d.ds=(rand()%100)/100.0,ptr[cnt++]=d;

#if 0
  exaSort(arr,exaScalar_t,offsetof(Data,ds),exaSortAlgoBinSort,1,
      exaGetComm(h));
  assert(exaArrayGetSize(arr)==N);
#endif

  array_free(&arr);

  comm_free(&c);
  MPI_Finalize();

  return 0;
}
