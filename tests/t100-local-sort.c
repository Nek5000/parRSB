#include <genmap-impl.h>
#include <sort-impl.h>
#include <time.h>

typedef struct{
  double ds;
  slong dl;
  sint proc;
} Data;

#define N 10

int main(int argc,char *argv[]){
  struct array arr; array_init(Data,&arr,N);

  srand(time(0));

  int i,cnt; Data d; Data *ptr=arr.ptr;
  for(i=cnt=0; i<N/2; i++){
      d.dl=rand()%100,d.ds=(rand()%100)/100.0,ptr[cnt++]=d;
      d.dl=rand()%100,ptr[cnt++]=d;
  }

#if 0
  exaSortArray2(arr,exaScalar_t,offsetof(Data,ds),
    exaLong_t,offsetof(Data,dl));
#endif

  ptr=arr.ptr;
  for(i=0; i<N-1; i++){
    assert(ptr[i].ds<=ptr[i+1].ds && "Field ds is not sorted");
    if(i%2==0)
      assert(ptr[i].dl<=ptr[i+1].dl && "Field dl is not sorted");
  }

  array_free(&arr);
  return 0;
}
