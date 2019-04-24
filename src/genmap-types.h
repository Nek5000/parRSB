//
// Header for MPI
//
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _GENMAP_TYPES_H_
#define _GENMAP_TYPES_H_
//
// Genmap int, long and scalar types
//
typedef long long GenmapLong;
typedef unsigned long long GenmapULong;
#define GenmapLongFormat "%lld"
#define GenmapULongFormat "%llu"
#define GENMAP_LONG MPI_LONG_LONG
#define GENMAP_UNSIGNED_LONG MPI_UNSIGNED_LONG_LONG

typedef int GenmapInt;
typedef unsigned int GenmapUInt;
#define GenmapIntFormat "%d"
#define GenmapUIntFormat "%u"
#define GENMAP_INT MPI_INT
#define GENMAP_UNSIGNED_INT MPI_UNSIGNED_INT

typedef double GenmapScalar;
#define GenmapScalarFormat "%lf"
#define GENMAP_SCALAR MPI_DOUBLE
//
// GenmapCommExternal
//
typedef MPI_Datatype GenmapDataType;
typedef MPI_Comm GenmapCommExternal;
//
// Genmap Pointer types
//
typedef struct GenmapComm_private *GenmapComm;
typedef struct GenmapHandle_private *GenmapHandle;
typedef struct GenmapVector_private *GenmapVector;
typedef struct GenmapElement_private *GenmapElements;
#if defined(PARRSB_GPU)
typedef struct parRSBKrylov_private* parRSBKrylov;
#endif

#endif

#ifdef __cplusplus
}
#endif
