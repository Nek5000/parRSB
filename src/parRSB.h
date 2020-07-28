#ifndef _PARRSB_H_
#define _PARRSB_H_

#define fparRSB_partMesh \
  FORTRAN_UNPREFIXED(fparrsb_partmesh,FPARRSB_PARTMESH)
void fparRSB_partMesh(int *part,long long *vtx,int *nel,
  int *nve,int *options,int *comm,int *err);

int parRSB_partMesh(int *part,long long *vtx,int nel,int nve,
  int *options,MPI_Comm comm);

#define MAXDIM 3
typedef struct {
  int proc;
  int orig;
  int seq;
  unsigned long long id;
  double coord[MAXDIM];
} elm_rcb;

#define fparRCB_partMesh \
  FORTRAN_UNPREFIXED(fparrcb_partmesh,FPARRCB_PARTMESH)
void fparRCB_partMesh(int *part,int *seq,double *vtx,int *nel,int *nv,
  int *options,int *comm,int *err);

int parRCB_partMesh  (int *part,int *seq,double *vtx,int  nel,int  nv,
  int *options,MPI_Comm comm);

int parRCB(struct comm *ci,struct array *a,int ndim);
void rcb_local(struct array *a,uint start,uint end,int ndim,buffer *buf);
#endif
