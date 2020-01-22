#ifndef _PARRSB_H_
#define _PARRSB_H_

#include "genmap-impl.h"
#include "exa-impl.h"

#define fparRSB_partMesh\
  FORTRAN_UNPREFIXED(fparrsb_partmesh,FPARRSB_PARTMESH)

#define MAXNV 8 /* maximum number of vertices per element */
typedef struct {
  int proc;
  long long id;
  int part;
  long long vtx[MAXNV];
} elm_data;

void fparRSB_partMesh(int *part, long long *vtx, int *nel,
                      int *nve, int *options, int *comm, int *err);

int parRSB_partMesh(int *part, long long *vtx, int nel, int nve,
  int *options,MPI_Comm comm);

#endif
