#ifndef _PARRSB_H_
#define _PARRSB_H_

#include "genmap-impl.h"
#include "gencon-impl.h"
#include "exa-impl.h"
#include "exa-memory.h"

#define MAXNV 8 /* maximum number of vertices per element */
typedef struct {
  int proc;
  long long id;
  int part;
  long long vtx[MAXNV];
} elm_data;

#define fparRSB_partMesh\
  FORTRAN_UNPREFIXED(fparrsb_partmesh,FPARRSB_PARTMESH)
void fparRSB_partMesh(int *part, long long *vtx, int *nel,
                      int *nve, int *options, int *comm,
                      int *err);

int parRSB_partMesh(int *part, long long *vtx, int nel, int nve,
  int *options,MPI_Comm comm);

#define fparRSB_findConnectivity\
  FORTRAN_UNPREFIXED(fparrsb_findconnectivity,\
    FPARRSB_FINDCONNECTIVITY)
void fparRSB_findConnectivity(double *coord,int *nel,int *nDim,
  long long *periodicInfo,int *nPeriodicFaces,long long *vertexId,
  double *tol,MPI_Fint *fcomm,int *err);

int parRSB_findConnectivity(double *coord,int nel,int nDim,
  long long *periodicInfo,int nPeriodicFaces,long long *vertexId,
  double tol,MPI_Comm comm);

#endif
