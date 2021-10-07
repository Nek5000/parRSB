#include <limits.h>
#include <math.h>

#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

#define MM 500

int project(genmap_vector x, struct laplacian *gl, mgData d, genmap_vector ri,
            int max_iter, struct comm *gsc, buffer *buf) {
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
  double *P = tcalloc(double, lelt *MM);
  double *W = tcalloc(double, lelt *MM);

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

// Input z should be orthogonal to 1-vector, have unit norm.
// RQI should not change z.
static int rqi(genmap_vector y, struct laplacian *gl, mgData d, genmap_vector z,
               struct comm *gsc, int max_iter, int grammian, buffer *bff) {
  assert(z->size == y->size);
  assert(z->size == gl->lelt);

  uint lelt = gl->lelt;

  // Grammian
  GenmapScalar *Z, *GZ, *M, *rhs, *v, *buf;
  GenmapMalloc(max_iter * lelt, &Z);
  GenmapMalloc(lelt, &GZ);
  GenmapMalloc(max_iter * max_iter, &M);
  GenmapMalloc(max_iter, &rhs);
  GenmapMalloc(max_iter, &v);
  GenmapMalloc(max_iter * max_iter, &buf);

  metric_tic(gsc, PROJECT);
  int ppfi = project(y, gl, d, z, 100, gsc, bff);
  metric_toc(gsc, PROJECT);

  slong out[2][1], bfr[2][1];
  slong in = lelt;
  comm_scan(out, gsc, gs_long, gs_add, &in, 1, bfr);
  slong nelg = out[1][0];

  genmap_vector err;
  genmap_vector_create(&err, lelt);

  uint i, j, k, l;
  for (i = 0; i < max_iter; i++) {
    GenmapScalar norm = genmap_vector_dot(y, y);
    comm_allreduce(gsc, gs_double, gs_add, &norm, 1, bfr);
    GenmapScalar normi = 1.0 / sqrt(norm);

    genmap_vector_axpby(z, z, 0.0, y, normi);
    genmap_vector_ortho_one(gsc, z, nelg);

    int N = i + 1;

    if (grammian == 1) {
      // if k>1;
      //  Z(:,k)=z-Z(:,1:k-1)*(Z(:,1:k-1)'*z);
      //  Z(:,k)=Z(:,k)/norm(Z(:,k));
      // end;
      if (i > 0) {
        // rhs = Z[1:k-1,:]*z
        for (j = 0; j < i; j++) {
          rhs[j] = 0.0;
          for (l = 0; l < lelt; l++)
            rhs[j] += Z[j * lelt + l] * z->data[l];
        }
        // Global reduction rhs[j]
        comm_allreduce(gsc, gs_double, gs_add, rhs, i, bfr);

        // Z[k,:] = z[:] - Z[:,1:lelt]*rhs[:]
        for (l = 0; l < lelt; l++)
          Z[i * lelt + l] = z->data[l];

        for (j = 0; j < i; j++) {
          for (l = 0; l < lelt; l++)
            Z[i * lelt + l] = Z[i * lelt + l] - rhs[j] * Z[j * lelt + l];
        }

        // Z[k,:]= Z[k,:]/||Z[k,:]||
        norm = 0.0;
        for (l = 0; l < lelt; l++)
          norm += Z[i * lelt + l] * Z[i * lelt + l];

        comm_allreduce(gsc, gs_double, gs_add, &norm, 1, bfr);
        norm = 1.0 / sqrt(norm);

        for (l = 0; l < lelt; l++)
          Z[i * lelt + l] *= norm;

        // M=Z(1:k,:)*G*Z(1:k,:);
        for (j = 0; j < N; j++) {
          GenmapLaplacianWeighted(GZ, gl, &Z[j * lelt], bff);
          for (k = 0; k < N; k++) {
            M[k * N + j] = 0.0;
            for (l = 0; l < lelt; l++)
              M[k * N + j] += Z[k * lelt + l] * GZ[l];
          }
        }

        // Global reduction of M
        comm_allreduce(gsc, gs_double, gs_add, M, N * N, buf);

        // Inverse power iterarion on M
        genmap_inverse_power(v, N, M, 0);

        for (j = 0; j < lelt; j++)
          z->data[j] = 0.0;

        for (j = 0; j < N; j++) {
          for (k = 0; k < lelt; k++)
            z->data[k] += Z[j * lelt + k] * v[j];
        }
        genmap_vector_ortho_one(gsc, z, nelg);
      } else {
        // Z(k,:) = z;
        for (l = 0; l < lelt; l++)
          Z[i * lelt + l] = z->data[l];
      }
    }

    metric_tic(gsc, PROJECT);
    ppfi = project(y, gl, d, z, 100, gsc, bff);
    metric_toc(gsc, PROJECT);

    genmap_vector_ortho_one(gsc, y, nelg);

    GenmapScalar lambda = genmap_vector_dot(y, z);
    comm_allreduce(gsc, gs_double, gs_add, &lambda, 1, bfr);

    genmap_vector_axpby(err, y, 1.0, z, -lambda);
    GenmapScalar norme = genmap_vector_dot(err, err);
    comm_allreduce(gsc, gs_double, gs_add, &norme, 1, bfr);
    norme = sqrt(norme);

    if (ppfi == 1)
      break;
  }

  GenmapFree(Z);
  GenmapFree(GZ);
  GenmapFree(M);
  GenmapFree(rhs);
  GenmapFree(v);
  GenmapFree(buf);

  genmap_destroy_vector(err);

  return i + 1;
}

static double sign(GenmapScalar a, GenmapScalar b) {
  GenmapScalar m = b >= 0.0 ? 1.0 : -1.0;
  return fabs(a) * m;
}

static int tqli(genmap_vector **eVectors, genmap_vector *eValues,
                genmap_vector diagonal, genmap_vector upper, int id) {
  assert(diagonal->size == upper->size + 1);

  GenmapInt n = diagonal->size;

  genmap_vector d, e;
  genmap_vector_create(&d, n);
  genmap_vector_create(&e, n);

  genmap_vector_copy(d, diagonal);
  genmap_vector_copy(e, upper);
  e->data[n - 1] = 0.0;

  /* Create the vector to store eigenvalues */
  genmap_vector_create(eValues, n);
  /* Init to identity */
  GenmapMalloc(n, eVectors);
  GenmapInt i;
  for (i = 0; i < n; i++) {
    genmap_vector_create_zeros(&(*eVectors)[i], n);
    (*eVectors)[i]->data[i] = 1.0;
  }

  GenmapInt j, k, l, iter, m;
  for (l = 0; l < n; l++) {
    iter = 0;
    do {
      for (m = l; m < n - 1; m++) {
        GenmapScalar dd = fabs(d->data[m]) + fabs(d->data[m + 1]);
        /* Should use a tolerance for this check */
        if (fabs(e->data[m]) / dd < GENMAP_TOL)
          break;
      }

      if (m != l) {
        if (iter++ == 30) {
          if (id == 0)
            printf("Too many iterations.\n");
          genmap_vector_copy(*eValues, d);
          return 1;
        }

        GenmapScalar g = (d->data[l + 1] - d->data[l]) / (2.0 * e->data[l]);
        GenmapScalar r = sqrt(g * g + 1.0);

        g = d->data[m] - d->data[l] + e->data[l] / (g + sign(r, g));
        GenmapScalar s = 1.0, c = 1.0, p = 0.0;

        for (i = m - 1; i >= l; i--) {
          GenmapScalar f = s * e->data[i];
          GenmapScalar b = c * e->data[i];

          if (fabs(f) >= fabs(g)) {
            c = g / f;
            r = sqrt(c * c + 1.0);
            e->data[i + 1] = f * r;
            s = 1.0 / r;
            c = c * s;
          } else {
            s = f / g;
            r = sqrt(s * s + 1.0);
            e->data[i + 1] = g * r;
            c = 1.0 / r;
            s = s * c;
          }

          g = d->data[i + 1] - p;
          r = (d->data[i] - g) * s + 2.0 * c * b;
          p = s * r;
          d->data[i + 1] = g + p;
          g = c * r - b;
          /* Find eigenvectors */
          for (k = 0; k < n; k++) {
            f = (*eVectors)[k]->data[i + 1];
            (*eVectors)[k]->data[i + 1] = s * (*eVectors)[k]->data[i] + c * f;
            (*eVectors)[k]->data[i] = c * (*eVectors)[k]->data[i] - s * f;
          }
          /* Done with eigenvectors */
        }

        if (r < GENMAP_TOL && i >= l)
          continue;

        d->data[l] -= p;
        e->data[l] = g;
        e->data[m] = 0.0;
      }
    } while (m != l);
  }

  /* Orthnormalize eigenvectors -- Just normalize? */
  for (i = 0; i < n; i++) {
    for (j = 0; j < i; j++) {
      GenmapScalar tmp = (*eVectors)[i]->data[j];
      (*eVectors)[i]->data[j] = (*eVectors)[j]->data[i];
      (*eVectors)[j]->data[i] = tmp;
    }
  }

  for (k = 0; k < n; k++) {
    e->data[k] = genmap_vector_dot((*eVectors)[k], (*eVectors)[k]);
    if (e->data[k] > 0.0)
      e->data[k] = sqrt(fabs(e->data[k]));
    GenmapScalar scale = 1.0 / e->data[k];
    genmap_vector_scale((*eVectors)[k], (*eVectors)[k], scale);
  }

  genmap_vector_copy(*eValues, d);

  genmap_destroy_vector(d);
  genmap_destroy_vector(e);

  return 0;
}

static int lanczos_aux(genmap_vector diag, genmap_vector upper,
                       genmap_vector **ri, uint lelt, int niter,
                       genmap_vector f, struct laplacian *gl, struct comm *gsc,
                       buffer *bfr) {
  assert(diag->size == niter);
  assert(diag->size == upper->size + 1);
  assert(f->size == lelt);

  slong out[2][1], buf[2][1];
  slong in = lelt;
  comm_scan(out, gsc, gs_long, gs_add, &in, 1, buf);
  slong start = out[0][0];
  slong nelg = out[1][0];

  if (nelg < niter) {
    niter = nelg;
    diag->size = niter;
    upper->size = niter - 1;
  }

  GenmapMalloc(niter + 1, ri);
  genmap_vector *rr = *ri;
  GenmapInt i;
  for (i = 0; i < niter + 1; ++i)
    rr[i] = NULL;

  genmap_vector r, p, w;
  genmap_vector_create_zeros(&p, lelt);
  genmap_vector_create(&w, lelt);
  genmap_vector_create(&r, lelt);

  genmap_vector_copy(r, f);

  genmap_vector_ortho_one(gsc, r, nelg);
  GenmapScalar rtr = genmap_vector_dot(r, r);
  comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, buf);
  GenmapScalar rnorm = sqrt(rtr);
  GenmapScalar rni = 1.0 / rnorm;

  genmap_vector_create(&rr[0], lelt);
  genmap_vector_scale(rr[0], r, rni);

  GenmapScalar eps = 1.e-5;
  GenmapScalar rtol = rnorm * eps;

  GenmapScalar rtz1 = 1.0;
  GenmapScalar pap = 0.0;
  GenmapScalar alpha, beta;
  GenmapScalar rtz2, pap_old;

  int iter;
  for (iter = 0; iter < niter; iter++) {
    rtz2 = rtz1;
    rtz1 = rtr;
    beta = rtz1 / rtz2;
    if (iter == 0)
      beta = 0.0;

    /* add2s1(p,r,beta,n) */
    for (i = 0; i < lelt; i++)
      p->data[i] = beta * p->data[i] + r->data[i];

    GenmapScalar pp = genmap_vector_dot(p, p);
    comm_allreduce(gsc, gs_double, gs_add, &pp, 1, buf);

    genmap_vector_ortho_one(gsc, p, nelg);

    metric_tic(gsc, LAPLACIAN);
    GenmapLaplacianWeighted(w->data, gl, p->data, bfr);
    metric_tic(gsc, LAPLACIAN);

    GenmapScalar ww = genmap_vector_dot(w, w);
    comm_allreduce(gsc, gs_double, gs_add, &ww, 1, buf);

    pap_old = pap;
    pap = genmap_vector_dot(w, p);
    comm_allreduce(gsc, gs_double, gs_add, &pap, 1, buf);
    // if (gsc->id == 0)
    //   printf("iter = %d beta = %lf pp = %lf pap = %lf ww = %lf\n", iter,
    //          beta, pp, pap, ww);

    alpha = rtz1 / pap;
    genmap_vector_axpby(r, r, 1.0, w, -1.0 * alpha);

    rtr = genmap_vector_dot(r, r);
    comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, buf);
    rnorm = sqrt(rtr);
    rni = 1.0 / rnorm;

    genmap_vector_create(&rr[iter + 1], lelt);
    genmap_vector_scale(rr[iter + 1], r, rni);

    if (iter == 0) {
      diag->data[iter] = pap / rtz1;
    } else {
      diag->data[iter] = (beta * beta * pap_old + pap) / rtz1;
      upper->data[iter - 1] = -beta * pap_old / sqrt(rtz2 * rtz1);
    }

    if (rnorm < rtol) {
      diag->size = iter + 1;
      upper->size = iter;
      iter = iter + 1;
      break;
    }
  }

  metric_acc(TOL_FINAL, rnorm);
  metric_acc(TOL_TARGET, rtol);

  genmap_destroy_vector(p);
  genmap_destroy_vector(w);
  genmap_destroy_vector(r);

  return iter;
}

static int lanczos(genmap_vector fiedler, uint lelt, int max_iter,
                   genmap_vector initv, struct laplacian *gl, struct comm *gsc,
                   buffer *buf) {
  genmap_vector alphaVec, betaVec;
  genmap_vector_create(&alphaVec, max_iter);
  genmap_vector_create(&betaVec, max_iter - 1);
  genmap_vector *q = NULL;

  metric_tic(gsc, LANCZOS);
  int iter =
      lanczos_aux(alphaVec, betaVec, &q, lelt, max_iter, initv, gl, gsc, buf);
  metric_toc(gsc, LANCZOS);

  genmap_vector evTriDiag;
  genmap_vector_create(&evTriDiag, iter);

  /* Use TQLI and find the minimum eigenvalue and associated vector */
  genmap_vector *eVectors, eValues;
  tqli(&eVectors, &eValues, alphaVec, betaVec, gsc->id);

  GenmapScalar eValMin = fabs(eValues->data[0]);
  GenmapInt eValMinI = 0;
  uint i;
  for (i = 1; i < iter; i++) {
    if (fabs(eValues->data[i]) < eValMin) {
      eValMin = fabs(eValues->data[i]);
      eValMinI = i;
    }
  }
  genmap_vector_copy(evTriDiag, eVectors[eValMinI]);

  GenmapInt j;
  for (i = 0; i < lelt; i++) {
    fiedler->data[i] = 0.0;
    for (j = 0; j < iter; j++)
      fiedler->data[i] += q[j]->data[i] * evTriDiag->data[j];
  }

  genmap_destroy_vector(alphaVec);
  genmap_destroy_vector(betaVec);
  genmap_destroy_vector(evTriDiag);
  genmap_destroy_vector(eValues);
  for (i = 0; i < iter; i++)
    genmap_destroy_vector(eVectors[i]);
  GenmapFree(eVectors);

  for (i = 0; i < iter + 1; i++)
    genmap_destroy_vector(q[i]);
  GenmapFree(q);

  return iter;
}

int fiedler(struct rsb_element *elems, uint lelt, int nv, int max_iter,
            int max_pass, int algo, struct comm *gsc, buffer *buf) {
  slong out[2][1], bfr[2][1];
  slong in = lelt;
  comm_scan(out, gsc, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0];
  slong nelg = out[1][0];

  struct csr_mat M;
  mgData d;
  if (algo == 1) {
    GenmapInitLaplacian(&M, elems, lelt, nv, gsc, buf);
    mgSetup(gsc, &M, &d);
  }

  struct laplacian gl;
  GenmapInitLaplacianWeighted(&gl, elems, lelt, nv, gsc, buf);

  genmap_vector initv, fiedler;
  genmap_vector_create(&initv, lelt);
  genmap_vector_create(&fiedler, lelt);

  uint i;
  for (i = 0; i < lelt; i++) {
    initv->data[i] = start + i + 1.0;
    if (start + i < nelg / 2)
      initv->data[i] += 1000 * nelg;
  }

  int ipass = 0, iter = 0, global = 1;
  do {
    genmap_vector_ortho_one(gsc, initv, nelg);

    double rtr = genmap_vector_dot(initv, initv);
    double rni;
    comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, &rni);
    rni = 1.0 / sqrt(rtr);
    genmap_vector_scale(initv, initv, rni);

    if (algo == 0)
      iter = lanczos(fiedler, lelt, max_iter, initv, &gl, gsc, buf);
    else if (algo == 1)
      iter = rqi(fiedler, &gl, d, initv, gsc, max_iter, 0 /* grammian */, buf);
    metric_acc(FIEDLER_NITER, iter);

    GenmapScalar norm = 0;
    for (i = 0; i < lelt; i++)
      norm += fiedler->data[i] * fiedler->data[i];
    GenmapScalar normi;
    comm_allreduce(gsc, gs_double, gs_add, &norm, 1, &normi);
    normi = 1.0 / sqrt(norm);

    genmap_vector_scale(fiedler, fiedler, normi);
    for (i = 0; i < lelt; i++)
      initv->data[i] = fiedler->data[i];
  } while (++ipass < max_pass && iter == max_iter);

  for (i = 0; i < lelt; i++)
    elems[i].fiedler = initv->data[i];

  genmap_destroy_vector(initv);
  genmap_destroy_vector(fiedler);
  laplacian_free(&gl);

  if (algo == 1) {
    mgFree(d);
    csr_mat_free(&M);
  }

  return 0;
}
