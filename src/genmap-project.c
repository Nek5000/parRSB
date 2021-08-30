#include <genmap-impl.h>

#define MM 500

int project(genmap_vector x, struct gs_laplacian *gl, mgData d,
            genmap_vector ri, int max_iter, struct comm *gsc, buffer *buf) {
  assert(x->size == ri->size);
  uint lelt = x->size;

  genmap_vector z0, z, dz, w, p, r;
  genmap_vector_create(&z, lelt);
  genmap_vector_create(&w, lelt);
  genmap_vector_create(&r, lelt);
  genmap_vector_create(&p, lelt);
  genmap_vector_create(&z0, lelt);
  genmap_vector_create(&dz, lelt);

  assert(max_iter < MM);
  double *P = tcalloc(double, lelt * MM);
  double *W = tcalloc(double, lelt * MM);

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
  genmap_vector_copy(p, z);

  GenmapScalar rz1 = genmap_vector_dot(r, z);
  comm_allreduce(gsc, gs_double, gs_add, &rz1, 1, bfr);

  GenmapScalar rr = genmap_vector_dot(r, r);
  comm_allreduce(gsc, gs_double, gs_add, &rr, 1, bfr);

  GenmapScalar alpha, beta, rz0, rz2, scale;

  double tol = 1e-5;
  double res_tol = rr * tol;

  uint j, k;
  i = 0;
  while (i < max_iter) {
    metric_tic(gsc, LAPLACIAN);
    GenmapLaplacianWeighted(w->data, gl, p->data, buf);
    metric_toc(gsc, LAPLACIAN);

    GenmapScalar den = genmap_vector_dot(p, w);
    comm_allreduce(gsc, gs_double, gs_add, &den, 1, bfr);
    alpha = rz1 / den;

    scale = 1.0 / sqrt(den);
    for (j = 0; j < lelt; j++) {
      W[i * lelt + j] = scale * w->data[j];
      P[i * lelt + j] = scale * p->data[j];
    }

    genmap_vector_axpby(x, x, 1.0, p, alpha);
    genmap_vector_axpby(r, r, 1.0, w, -alpha);

    rr = genmap_vector_dot(r, r);
    comm_allreduce(gsc, gs_double, gs_add, &rr, 1, bfr);

    if (rr < res_tol || sqrt(rr) < tol)
      break;

    GenmapScalar norm0 = genmap_vector_dot(z, z);
    comm_allreduce(gsc, gs_double, gs_add, &norm0, 1, bfr);

    genmap_vector_copy(z0, z);
    mg_vcycle(z->data, r->data, d);

    GenmapScalar norm1 = genmap_vector_dot(z, z);
    comm_allreduce(gsc, gs_double, gs_add, &norm1, 1, bfr);

    rz0 = rz1;
    genmap_vector_ortho_one(gsc, z, nelg);
    rz1 = genmap_vector_dot(r, z);
    comm_allreduce(gsc, gs_double, gs_add, &rz1, 1, bfr);

    genmap_vector_axpby(dz, z, 1.0, z0, -1.0);
    rz2 = genmap_vector_dot(r, dz);
    comm_allreduce(gsc, gs_double, gs_add, &rz2, 1, bfr);

    beta = rz2 / rz0;
    genmap_vector_axpby(p, z, 1.0, p, beta);

    i++;

    metric_tic(gsc, PROJECT);
    for (k = 0; k < lelt; k++)
      P[(MM - 1) * lelt + k] = 0.0;

    for (j = 0; j < i; j++) {
      double a = 0.0;
      for (k = 0; k < lelt; k++)
        a += W[j * lelt + k] * p->data[k];
      comm_allreduce(gsc, gs_double, gs_add, &a, 1, bfr);
      for (k = 0; k < lelt; k++)
        P[(MM - 1) * lelt + k] += a * P[j * lelt + k];
    }

    for (k = 0; k < lelt; k++)
      p->data[k] -= P[(MM - 1) * lelt + k];
    metric_toc(gsc, PROJECT);
  }

  genmap_destroy_vector(z);
  genmap_destroy_vector(w);
  genmap_destroy_vector(p);
  genmap_destroy_vector(r);
  genmap_destroy_vector(z0);
  genmap_destroy_vector(dz);

  GenmapFree(P);
  GenmapFree(W);

  return i + 1;
}
