#ifndef _PARRSB_SORT_H_
#define _PARRSB_SORT_H_

#include <gslib.h>

#if !defined(MPI)
#error "gslib needs to be compiled with MPI"
#endif

#if !defined(GLOBAL_LONG_LONG)
#error "gslib needs to be compiled with GLOBAL_LONG_LONG"
#endif

#ifdef __cplusplus
extern "C" {
#endif

int parallel_sort_(struct array *array, unsigned algo, unsigned balance,
                   size_t usize, size_t align, struct comm *c, buffer *bfr,
                   unsigned nfields, ...);

#define parallel_sort(T, array, field, type, algo, balance, c, bfr)            \
  parallel_sort_(array, algo, balance, sizeof(T), ALIGNOF(T), c, bfr, 1,       \
                 offsetof(T, field), type)

#define parallel_sort_2(T, array, f1, t1, f2, t2, algo, balance, c, bfr)       \
  parallel_sort_(array, algo, balance, sizeof(T), ALIGNOF(T), c, bfr, 2,       \
                 offsetof(T, f1), t1, offsetof(T, f2), t2)

#ifdef __cplusplus
}
#endif

#endif
