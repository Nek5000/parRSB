#ifndef _SORT_IMPL_H_
#define _SORT_IMPL_H_

#include <genmap-impl.h>

#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)>(b) ? (a) : (b))
//
// exaSort: general functions
//
typedef enum{
  exaSortAlgoBinSort      =0,
  exaSortAlgoHyperCubeSort=1
} exaSortAlgo;

typedef struct{
  int nfields;
  gs_dom t[3];
  uint offset[3];

  struct array *a;
  size_t unit_size,align;

  int balance;
  exaSortAlgo algo;

  buffer buf;
} sort_data_private;
typedef sort_data_private* sort_data;

double get_scalar(struct array *a,uint i,uint offset,uint usize,
  gs_dom type);
void get_extrema(void *extrema_,sort_data data,uint field,struct comm *c);
int  set_dest(uint *proc,uint np,ulong start,uint size,ulong nelem);
int  load_balance(struct array *a,size_t size,struct comm *c,
  struct crystal *cr);
int  sort_local(sort_data data);
int  sort_private(sort_data data,struct comm *c);
//
// exaBinSort
//
int exaBinSort(sort_data data,struct comm *c);
//
// exaHyperCubeSort
//
typedef struct{
  sort_data data;
  int nProbes;
  double *probes;
  slong *probe_cnt;
} hypercube_sort_data_private;
typedef hypercube_sort_data_private* hypercube_sort_data;

int exaHyperCubeSort(hypercube_sort_data data,struct comm *c);

#define parallel_sort(T,A,off,type,c) do {\
  sort_data_private sd;\
  sd.unit_size=sizeof(T);\
  sd.align=ALIGNOF(T);\
  sd.nfields=1;\
  sd.t[0]=type;\
  sd.offset[0]=off;\
  sd.a=A;\
  sd.algo=exaSortAlgoBinSort;\
  sd.balance=1;\
  buffer_init(&sd.buf,1024);\
  sort_private(&sd,c);\
  buffer_free(&sd.buf);\
} while (0)

#endif // exasort-impl
