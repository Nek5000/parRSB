#ifndef _GENMAP_GSLIB_H_
#define _GENMAP_GSLIB_H_

#ifdef __cplusplus
extern "C" {
#endif

// Data type sint/uint
// (defualt) int
// #define USE_LONG long
// #define USE_LONG_LONG long long

// Data type slong/ulong
// (default) int
// #define GLOBAL_LONG long
// #define GLOBAL_LONG_LONG long long

#include "gslib.h"

#if !defined(MPI)
#error "gslib needs to be compiled with MPI"
#endif

#if !defined(GLOBAL_LONG_LONG)
#error "gslib needs to be compiled with GLOBAL_LONG_LONG"
#endif

#define TYPE_INT    0
#define TYPE_LONG   1
#define TYPE_FLOAT  2
#define TYPE_DOUBLE 3

#define genmap_gs_long gs_long_long
#define genmap_gs_int gs_int
#define genmap_gs_scalar gs_double

#endif

#ifdef __cplusplus
}
#endif
