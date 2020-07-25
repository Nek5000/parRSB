#ifndef _GENCON_H_
#define _GENCON_H_

#include "exa.h"
#include "exasort.h"

/* Upper bound for number of dimensions */
#define GC_MAX_DIM 3

/* Boundary condition types */
#define GC_PERIODIC "P  "

/* Upper bounds for elements */
#define GC_MAX_FACES 6
#define GC_MAX_VERTICES 8
#define GC_MAX_NEIGHBORS 3

/* Upper bounds for faces */
#define GC_MAX_FACE_VERTICES 4

/* Header lengths */
#define GC_RE2_HEADER_LEN 80
#define GC_CO2_HEADER_LEN 132

typedef struct Point_private *Point;
typedef struct Element_private *Element;
typedef struct Boundary_private *BoundaryFace;
typedef struct Mesh_private *Mesh;

extern int NEIGHBOR_MAP[GC_MAX_VERTICES][GC_MAX_NEIGHBORS];
extern int PRE_TO_SYM_VERTEX[GC_MAX_VERTICES];
extern int PRE_TO_SYM_FACE[GC_MAX_FACES];

/* Mesh */
int MeshInit(Mesh *m_,int nel,int nDim);
Element MeshGetElements(Mesh m);
int MeshFree(Mesh m);

int read_re2_mesh(Mesh *mesh,char *fname,struct comm *c);
int read_co2_file(Mesh  mesh,char *fname,struct comm *c);
int write_co2_file(Mesh mesh,char *fname,struct comm *c);
int read_co2_mesh(Mesh *mesh,char *fname,struct comm *c);

int findMinNeighborDistance(exaHandle h,Mesh mesh);
int findSegments(exaHandle h,Mesh mesh,exaScalar tol);
int setGlobalID(exaHandle h,Mesh mesh);
int sendBack(exaHandle h,Mesh mesh);
int matchPeriodicFaces(exaHandle h,Mesh mesh);

#endif
