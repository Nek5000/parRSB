#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <gencon-impl.h>
#include <parRSB.h>

void fparRSB_findConnectivity(long long *vertexId,double *coord,
  int *nelt,int *ndim,long long *periodicInfo,int *nPeriodicFaces,
  double *tol,MPI_Fint *fcomm,int *err)
{
  *err=1;
  MPI_Comm c=MPI_Comm_f2c(*fcomm);
  *err=parRSB_findConnectivity(vertexId,coord,*nelt,*ndim,periodicInfo,
      *nPeriodicFaces,*tol,c);
}

int transferBoundaryFaces_(Mesh mesh,struct comm *c){
  uint size=c->np;

  struct array *boundary=&mesh->boundary;
  BoundaryFace ptr=boundary->ptr;
  int nFaces=boundary->n;

  slong nelgt=mesh->nelgt;
  sint nelt=nelgt/size,nrem=nelgt-nelt*size;
  slong N=(size-nrem)*nelt;

  sint i; slong eid;
  for(i=0;i<nFaces;i++,ptr++){
    eid=ptr->elementId;
    if(eid<N) ptr->proc=eid/nelt;
    else ptr->proc=(eid-N)/(nelt+1)+size-nrem;
  }

  struct crystal cr; crystal_init(&cr,c);
  sarray_transfer(struct Boundary_private,boundary,proc,1,&cr);
  crystal_free(&cr);
}

// coord[nelt,nv,ndim] - in, vertices are orders in preprocessor ordering
// vertexid[nelt,nv] - out
int parRSB_findConnectivity(long long *vertexid,double *coord,
  int nelt,int ndim,long long *periodicInfo,int nPeriodicFaces,
  double tol,MPI_Comm comm)
{
  struct comm c; comm_init(&c,comm);
  uint rank=c.id,size=c.np;

  Mesh mesh; MeshInit(&mesh,nelt,ndim);

  slong out[2][1],buff[2][1],in=nelt;
  comm_scan(out,&c,gs_long,gs_add,&in,1,buff);
  ulong start=out[0][0];
  ulong nelgt=mesh->nelgt=mesh->nelgv=out[1][0];

  int nvertex=mesh->nVertex;
  uint nunits=nvertex*nelt;

  struct Point_private p;
  uint i,j,k,l;
  for(i=0; i<nelt; i++){
    for(k=0; k<nvertex; k++){
      j=PRE_TO_SYM_VERTEX[k];
      for(l=0; l<ndim; l++)
        p.x[l]=coord[i*nvertex*ndim+j*ndim+l];
      p.elementId =start+i;
      p.sequenceId=nvertex*(start+i)+k;
      p.origin    =rank;

      array_cat(struct Point_private,&mesh->elements,&p,1);
    }
  }
  assert(mesh->elements.n==nunits);

  struct Boundary_private b;
  for(i=0; i<nPeriodicFaces; i++){
    b.elementId=periodicInfo[4*i+0]-1;
    b.faceId   =PRE_TO_SYM_FACE[periodicInfo[4*i+1]-1];
    b.bc[0]    =periodicInfo[4*i+2]-1;
    b.bc[1]    =PRE_TO_SYM_FACE[periodicInfo[4*i+3]-1];
    array_cat(struct Boundary_private,&mesh->boundary,&b,1);
    //printf("eid/faceid/bc[0]/bc[1]:%lld %d %lld %d\n",
    //  b.elementId,b.faceId,b.bc[0],b.bc[1]);
  }
  assert(mesh->boundary.n==nPeriodicFaces);

  transferBoundaryFaces(mesh,&c);

  findMinNeighborDistance(mesh);
  findSegments(mesh,&c,tol);
  setGlobalID(mesh,&c);
  sendBack(mesh,&c);
  matchPeriodicFaces(mesh,&c);

  // copy output
  Point ptr=mesh->elements.ptr;
  k=0;
  for(i=0; i<nelt; i++){
    vertexid[k++]=ptr->elementId+1;
    for(j=0; j<nvertex; j++){
      vertexid[k++]=ptr->globalId+1;
      ptr++;
    }
  }

  MeshFree(mesh);
  comm_free(&c);

  return 0;
}
