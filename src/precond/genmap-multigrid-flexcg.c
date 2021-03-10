#include <math.h>
#include <stdio.h>

#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

int flex_cg(genmap_handle h, struct comm *gsc, mgData d, genmap_vector ri,
            int maxIter, genmap_vector x) {
  assert(x->size == ri->size);
  assert(x->size == genmap_get_nel(h));

  uint lelt = x->size;
  GenmapLong nelg = genmap_get_partition_nel(h);

  genmap_vector z0, z, dz, w, p, r;
  GenmapCreateVector(&z, lelt);
  GenmapCreateVector(&w, lelt);
  GenmapCreateVector(&r, lelt);
  GenmapCreateVector(&p, lelt);
  GenmapCreateVector(&z0, lelt);
  GenmapCreateVector(&dz, lelt);

#define PREC 1
#define ORTH 1

  int rank = gsc->id;

  uint i;
  for (i = 0; i < lelt; i++)
    x->data[i] = 0.0, r->data[i] = ri->data[i];

  GenmapCopyVector(z, r);
#if ORTH
  GenmapOrthogonalizebyOneVector(gsc, z, nelg);
#endif

  GenmapScalar den, alpha, beta, rz0, rz1, rz2, rr, buf;

  rz1 = GenmapDotVector(r, z);
  comm_allreduce(gsc, gs_double, gs_add, &rz1, 1, &buf);

  GenmapCopyVector(p, z);

  i = 0;
  while (i < maxIter && sqrt(rz1) > GENMAP_TOL) {
    GenmapLaplacian(h, p->data, w->data);

    den = GenmapDotVector(p, w);
    comm_allreduce(gsc, gs_double, gs_add, &den, 1, &buf);

    alpha = rz1 / den;

    GenmapAxpbyVector(x, x, 1.0, p, alpha);
    GenmapAxpbyVector(r, r, 1.0, w, -alpha);

    GenmapCopyVector(z0, z);
#if PREC
    mg_vcycle(z->data, r->data, d);
#else
    GenmapCopyVector(z, r);
#endif
#if ORTH
    GenmapOrthogonalizebyOneVector(gsc, z, nelg);
#endif

    rz0 = rz1;

    rz1 = GenmapDotVector(r, z);
    comm_allreduce(gsc, gs_double, gs_add, &rz1, 1, &buf);

    GenmapAxpbyVector(dz, z, 1.0, z0, -1.0);
    rz2 = GenmapDotVector(r, dz);
    comm_allreduce(gsc, gs_double, gs_add, &rz2, 1, &buf);
    beta = rz2 / rz0;

    GenmapAxpbyVector(p, z, 1.0, p, beta);
    i++;

    rr = GenmapDotVector(r, r);
    comm_allreduce(gsc, gs_double, gs_add, &rr, 1, &buf);
  }

  GenmapDestroyVector(z);
  GenmapDestroyVector(w);
  GenmapDestroyVector(p);
  GenmapDestroyVector(r);
  GenmapDestroyVector(z0);
  GenmapDestroyVector(dz);

  return i;
}
