#include <exasort.h>
#include <exa-memory.h>
#include <time.h>

typedef struct{
  exaScalar ds;
  slong dl;
  uint proc;
} Data;

#define N 10

void check(exaArray arr,exaHandle h){
  uint size=exaSize(h);
  uint rank=exaRank(h);

  Data *ptr=exaArrayGetPointer(arr);
  uint n=exaArrayGetSize(arr),i;
  for(i=0; i<n-1; i++){
    assert(ptr[i].ds<=ptr[i+1].ds && "Field ds is not sorted");
    if(i%2==0)
      assert(ptr[i].dl<=ptr[i+1].dl && "Field dl is not sorted");
  }

  uint root=0;

  struct array a; array_init(Data,&a,2); Data *p=a.ptr;
  p[0]     =ptr[0],p[1]     =ptr[n-1];
  p[0].proc=root  ,p[1].proc=root;
  p[0].dl  =2*rank,p[1].dl  =2*rank+1;
  a.n=2;

  struct comm c; comm_init(&c,MPI_COMM_WORLD);
  struct crystal cr; crystal_init(&cr,&c);

  sarray_transfer(Data,&a,proc,0,&cr);
  if(rank==root)
    assert(a.n==2*size);

  buffer buf; buffer_init(&buf,1024);
  sarray_sort(Data,a.ptr,a.n,dl,1,&buf);
  buffer_free(&buf);

  ptr=a.ptr;
  if(rank==root)
    for(i=0; i<a.n-1; i++)
      assert(ptr[i].ds<=ptr[i+1].ds && "Field ds is not sorted globally");

  crystal_free(&cr);
  comm_free(&c);
  array_free(&a);
}

int main(int argc,char *argv[]){
  MPI_Init(&argc,&argv);

  exaHandle h;
  exaInit(&h,MPI_COMM_WORLD,"/host");

  srand(time(0));

  exaArray arr; exaArrayInit(&arr,Data,N);

  int i; Data d;
  for(i=0; i<N/2; i++){
      d.dl=rand()%100,d.ds=(rand()%100)/100.0,exaArrayAppend(arr,&d);
      d.dl=rand()%100,exaArrayAppend(arr,&d);
  }

  exaSort2(arr,exaScalar_t,offsetof(Data,ds),
    exaLong_t,offsetof(Data,dl),
    exaSortAlgoBinSort,1,exaGetComm(h));
  check(arr,h);

  exaArraySetSize(arr,0);
  for(i=0; i<N/2; i++){
      d.dl=rand()%100,d.ds=(rand()%100)/100.0,exaArrayAppend(arr,&d);
      d.dl=rand()%100,exaArrayAppend(arr,&d);
  }

  exaSort2(arr,exaScalar_t,offsetof(Data,ds),
    exaLong_t,offsetof(Data,dl),
    exaSortAlgoHyperCubeSort,1,exaGetComm(h));
  check(arr,h);

  exaDestroy(arr);
  exaFinalize(h);
  MPI_Finalize();

  return 0;
}
