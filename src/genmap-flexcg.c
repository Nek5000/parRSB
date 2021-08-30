#include <genmap-impl.h>

int flex_cg(genmap_vector x, struct gs_laplacian *gl, mgData d,
            genmap_vector ri, int maxIter, struct comm *gsc, buffer *buf) {
  assert(x->size == ri->size);
  uint lelt = x->size;

  genmap_vector z0, z, dz, w, p, r;
  genmap_vector_create(&z, lelt);
  genmap_vector_create(&w, lelt);
  genmap_vector_create(&r, lelt);
  genmap_vector_create(&p, lelt);
  genmap_vector_create(&z0, lelt);
  genmap_vector_create(&dz, lelt);

#define PREC 1
#define ORTH 1

  uint i;
  for (i = 0; i < lelt; i++) {
    x->data[i] = 0.0;
    r->data[i] = ri->data[i];
  }

  slong out[2][1], bfr[2][1];
  slong in = lelt;
  comm_scan(out, gsc, gs_long, gs_add, &in, 1, bfr);
  slong nelg = out[1][0];

  genmap_vector_copy(z, r);
#if ORTH
  genmap_vector_ortho_one(gsc, z, nelg);
#endif

  GenmapScalar den, alpha, beta, rz0, rz1, rz2, rr;

  rz1 = genmap_vector_dot(r, z);
  comm_allreduce(gsc, gs_double, gs_add, &rz1, 1, bfr);

  genmap_vector_copy(p, z);

  i = 0;
  while (i < maxIter && sqrt(rz1) > GENMAP_TOL) {
    metric_tic(gsc, LAPLACIAN);
    GenmapLaplacianWeighted(w->data, gl, p->data, buf);
    metric_toc(gsc, LAPLACIAN);

    den = genmap_vector_dot(p, w);
    comm_allreduce(gsc, gs_double, gs_add, &den, 1, bfr);

    alpha = rz1 / den;

    genmap_vector_axpby(x, x, 1.0, p, alpha);
    genmap_vector_axpby(r, r, 1.0, w, -alpha);

    genmap_vector_copy(z0, z);
#if PREC
    mg_vcycle(z->data, r->data, d);
#else
    genmap_vector_copy(z, r);
#endif
#if ORTH
    genmap_vector_ortho_one(gsc, z, nelg);
#endif

    rz0 = rz1;

    rz1 = genmap_vector_dot(r, z);
    comm_allreduce(gsc, gs_double, gs_add, &rz1, 1, bfr);

    genmap_vector_axpby(dz, z, 1.0, z0, -1.0);
    rz2 = genmap_vector_dot(r, dz);
    comm_allreduce(gsc, gs_double, gs_add, &rz2, 1, bfr);
    beta = rz2 / rz0;

    genmap_vector_axpby(p, z, 1.0, p, beta);
    i++;

    rr = genmap_vector_dot(r, r);
    comm_allreduce(gsc, gs_double, gs_add, &rr, 1, bfr);
  }

  genmap_destroy_vector(z);
  genmap_destroy_vector(w);
  genmap_destroy_vector(p);
  genmap_destroy_vector(r);
  genmap_destroy_vector(z0);
  genmap_destroy_vector(dz);

  return i;
}
