#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

#include <math.h>
#include <stdio.h>

#if 0
//input r should have zero-sum
int rqi(GenmapHandle h,GenmapComm c,GenmapVector r,GenmapVector fiedler){
  assert(r->size==fiedler->size);

  GenmapVector y; GenmapCreateVector(&y,r->size);
  flex_cg(h,c,d,r,

  GenmapDestroyVector(y);

  return iter;
}
#endif
