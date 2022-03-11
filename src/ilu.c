#include "ilu.h"
#include <math.h>

struct dof {
  slong v;
  uint src, dest;
}; 

static int find_keys(uint nelt, int nv, const slong *vtx, struct comm *c,
                     struct crystal *cr, buffer *bfr) { 
  uint ndofs = nelt * nv;
  int active = 0;
  if (ndofs > 0)
    active = 1;

  struct array dofs;
  array_init(struct dof, &dofs, npts);

  struct dof t = {0, c->id, 0};
  for (uint i = 0; i < ndofs; i++) {
    for (uint j = 0; j < nv; j++) {
      t.v = vtx[i * nv + j], t.sest = t.v % c->np;
      array_cat(struct dof, &dofs, &t, 1);
    }
  }

  sarray_transfer(struct dof, &dofs, dest, 1, cr);

  array_free(&dofs);

  return 0;
}
