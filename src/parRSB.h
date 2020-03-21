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
} elm_rsb;

#define MAXDIM 3
typedef struct {
  int proc;
  long long id;
  int part;
  double coord[MAXDIM];
} elm_rcb;

#define fparRSB_partMesh\
  FORTRAN_UNPREFIXED(fparrsb_partmesh,FPARRSB_PARTMESH)
void fparRSB_partMesh(int *part, long long *vtx, int *nel,
                      int *nve, int *options, int *comm,
                      int *err);

int parRSB_partMesh(int *part, long long *vtx, int nel, int nve,
  int *options,MPI_Comm comm);

#define fparRCB_partMesh\
  FORTRAN_UNPREFIXED(fparrcb_partmesh,FPARRCB_PARTMESH)
void fparRCB_partMesh(int *part,double *vtx,int *nel, int *ndim,
  int *options,int *comm,int *err);

int parRCB_partMesh(int *part,double *vtx,int nel,int ndim,
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
