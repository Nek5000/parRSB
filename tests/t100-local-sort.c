#include <exasort.h>
#include <time.h>

typedef struct{
  exaScalar ds;
  slong dl;
  sint proc;
} Data;

#define N 10

int main(int argc,char *argv[]){
  exaArray arr;
  exaArrayInit(&arr,Data,1);

  srand(time(0));

  int i; Data d;
  for(i=0; i<N/2; i++){
      d.dl=rand()%100,d.ds=(rand()%100)/100.0,exaArrayAppend(arr,&d);
      d.dl=rand()%100,exaArrayAppend(arr,&d);
  }

  exaSortArray2(arr,exaScalar_t,offsetof(Data,ds),
    exaLong_t,offsetof(Data,dl));

  Data *ptr=exaArrayGetPointer(arr);
  for(i=0; i<N-1; i++){
    assert(ptr[i].ds<=ptr[i+1].ds && "Field ds is not sorted");
    if(i%2==0)
      assert(ptr[i].dl<=ptr[i+1].dl && "Field dl is not sorted");
  }

  exaDestroy(arr);

  return 0;
}
