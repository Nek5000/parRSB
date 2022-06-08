#include "coarse.h"
#include "genmap-impl.h"
#include "multigrid.h"
#include "sort.h"
#include <math.h>
#include <time.h>

#define MM 500

extern void matrix_inverse(int N, double *A);

inline static scalar dot(scalar *y, scalar *x, uint n) {
  scalar result = 0.0;
  for (uint i = 0; i < n; i++)
    result += x[i] * y[i];

  return result;
}

inline static void ortho(scalar *q, uint lelt, ulong n, struct comm *c) {
  uint i;
  scalar sum = 0.0;
  for (i = 0; i < lelt; i++)
    sum += q[i];

  scalar buf;
  comm_allreduce(c, gs_double, gs_add, &sum, 1, &buf);
  sum /= n;

  for (i = 0; i < lelt; i++)
    q[i] -= sum;
}

struct fiedler {
  GenmapScalar fiedler;
  uint proc, seq;
  int part0, part1;
};

static sint converged(struct array *fdlr, GenmapScalar *fiedler, struct comm *c,
                      slong nelg, buffer *buf) {
  struct fiedler *ptr = fdlr->ptr;

  uint i;
  for (i = 0; i < fdlr->n; i++)
    ptr[i].fiedler = fiedler[i];

  parallel_sort(struct fiedler, fdlr, fiedler, gs_double, 0, 1, c, buf);

  slong out[2][1], bfr[2][1], in = fdlr->n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0];

  ptr = fdlr->ptr;
  for (i = 0; i < fdlr->n; i++)
    if (start + i < nelg / 2)
      ptr[i].part1 = 0;
    else
      ptr[i].part1 = 1;

  sint changed = 0;
  for (i = 0; i < fdlr->n; i++) {
    if (ptr[i].part1 != ptr[i].part0)
      changed = 1;
    ptr[i].part0 = ptr[i].part1;
  }

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct fiedler, fdlr, proc, 0, &cr);
  sarray_sort(struct fiedler, fdlr->ptr, fdlr->n, seq, 0, buf);
  crystal_free(&cr);

  comm_allreduce(c, gs_int, gs_add, &changed, 1, bfr);
  return changed > 0 ? 0 : 1;
}

int power_serial(double *y, int N, double *A, int verbose) {
  time_t t;
  srand((unsigned)time(&t));

  int i;
  GenmapScalar norm = 0.0;
  for (i = 0; i < N; i++) {
    y[i] = (rand() % 50) / 50.0;
    norm += y[i] * y[i];
  }

  GenmapScalar normi = 1.0 / sqrt(norm);
  for (i = 0; i < N; i++)
    y[i] *= normi;

  double *Ay;
  GenmapCalloc(N, &Ay);

  int j, k, l;
  GenmapScalar err = 1.0, lambda;
  for (i = 0; i < 100; i++) {
    norm = 0.0;
    for (j = 0; j < N; j++) {
      Ay[j] = 0.0;
      for (k = 0; k < N; k++) {
        Ay[j] += A[j * N + k] * y[k];
      }
      norm += Ay[j] * Ay[j];
    }

    if (i > 0)
      err = (sqrt(norm) - lambda) / lambda;
    lambda = sqrt(norm);

    normi = 1.0 / sqrt(norm);
    for (j = 0; j < N; j++)
      y[j] = Ay[j] * normi;

    if (fabs(err) < 1.e-12)
      break;
  }

  GenmapFree(Ay);

  return i;
}

int inv_power_serial(double *y, int N, double *A, int verbose) {
  double *Ainv;
  GenmapCalloc(N * N, &Ainv);

  int j, k;
  for (j = 0; j < N; j++) {
    for (k = 0; k < N; k++)
      Ainv[j * N + k] = A[k * N + j];
  }

  matrix_inverse(N, Ainv);

  for (j = 0; j < N; j++) {
    for (k = 0; k < N; k++)
      A[j * N + k] = Ainv[k * N + j];
  }

  j = power_serial(y, N, Ainv, verbose);

  GenmapFree(Ainv);

  return j;
}

static int project(scalar *x, uint n, scalar *b, struct laplacian *L,
                   struct mg *d, struct comm *c, int miter, int null_space,
                   int verbose, buffer *bfr) {
  slong out[2][1], buf[2][1], in = n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
  ulong ng = out[1][0];

  if (ng == 0)
    return 0;

  scalar *z = (scalar *)tcalloc(scalar, 6 * n);
  scalar *w = z + n, *r = w + n, *p = r + n, *z0 = p + n, *dz = z0 + n;
  scalar *P = tcalloc(scalar, 2 * (miter + 1) * n);
  scalar *W = P + n * (miter + 1);

  uint i;
  for (i = 0; i < n; i++)
    x[i] = 0, r[i] = b[i];

  scalar rr = dot(r, r, n);
  comm_allreduce(c, gs_double, gs_add, &rr, 1, buf);
  scalar tol = 1e-3, rtol = rr * tol * tol;

  for (i = 0; i < n; i++)
    z[i] = r[i];
  if (null_space)
    ortho(z, n, ng, c);
  scalar rz1 = dot(z, z, n);
  comm_allreduce(c, gs_double, gs_add, &rz1, 1, buf);

  for (i = 0; i < n; i++)
    p[i] = z[i];

  scalar alpha, beta, rzt, rz2;

  uint j, k;
  for (i = 0; i < miter; i++) {
    // mat_vec_csr(w, p, S, gsh, wrk, bfr);
    laplacian(w, L, p, bfr);

    scalar pw = dot(p, w, n);
    comm_allreduce(c, gs_double, gs_add, &pw, 1, buf);
    alpha = rz1 / pw;

    pw = 1 / sqrt(pw);
    for (j = 0; j < n; j++)
      W[i * n + j] = pw * w[j], P[i * n + j] = pw * p[j];

    for (j = 0; j < n; j++)
      x[j] += alpha * p[j], r[j] -= alpha * w[j];

    rr = dot(r, r, n);
    comm_allreduce(c, gs_double, gs_add, &rr, 1, buf);

    if (rr < rtol || sqrt(rr) < tol)
      break;

    for (j = 0; j < n; j++)
      z0[j] = z[j];

#if 1
    mg_vcycle(z, r, d, c, bfr);
#else
    for (j = 0; j < n; j++)
      z[j] = r[j];
#endif

    rzt = rz1;
    if (null_space)
      ortho(z, n, ng, c);
    rz1 = dot(r, z, n);
    comm_allreduce(c, gs_double, gs_add, &rz1, 1, buf);

    for (j = 0; j < n; j++)
      dz[j] = z[j] - z0[j];
    rz2 = dot(r, dz, n);
    comm_allreduce(c, gs_double, gs_add, &rz2, 1, buf);

    if (c->id == 0 && verbose > 0)
      printf("rr = %lf rtol = %lf rz0 = %lf rz1 = %lf rz2 = %lf\n", rr, rtol,
             rzt, rz1, rz2);

    beta = rz2 / rzt;
    for (j = 0; j < n; j++)
      p[j] = z[j] + beta * p[j];

    for (k = 0; k < n; k++)
      P[miter * n + k] = 0;

    for (j = 0; j <= i; j++) {
      pw = 0;
      for (k = 0; k < n; k++)
        pw += W[j * n + k] * p[k];
      comm_allreduce(c, gs_double, gs_add, &pw, 1, buf);
      for (k = 0; k < n; k++)
        P[miter * n + k] += pw * P[j * n + k];
    }

    for (k = 0; k < n; k++)
      p[k] -= P[miter * n + k];
  }

  free(z);
  free(P);

  return i == miter ? i : i + 1;
}

// Input z should be orthogonal to 1-vector, have unit norm.
// inverse iteration should not change z.
static int inverse(scalar *y, uint lelt, struct rsb_element *elems, int nv,
                   scalar *z, struct comm *gsc, int miter, slong nelg,
                   buffer *buf) {
  struct laplacian *wl = laplacian_init(elems, lelt, nv, GS, gsc, buf);

  // Reserve enough memory in buffer
  size_t wrk = sizeof(ulong) * lelt + sizeof(slong) * nv * lelt;
  buffer_reserve(buf, wrk);

  slong out[2][1], bfr[2][1], in = lelt;
  comm_scan(out, gsc, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0];

  ulong *eid = (ulong *)buf->ptr;
  slong *vtx = (slong *)(eid + lelt);
  uint i, j, k, l;
  for (i = k = 0; i < lelt; i++) {
    eid[i] = start + i + 1;
    for (j = 0; j < nv; j++)
      vtx[k++] = elems[i].vertices[j];
  }

  struct crystal cr;
  crystal_init(&cr, gsc);
  struct par_mat *L = par_csr_setup_con(lelt, eid, vtx, nv, 1, gsc, &cr, buf);

  int factor = 2;
  char *fc = getenv("PARRSB_RSB_MG_FACTOR");
  if (fc != NULL)
    factor = atoi(fc);
  struct mg *d = mg_setup(L, factor, gsc, &cr, buf);
  crystal_free(&cr);

  scalar *err = tcalloc(scalar, 2 * lelt + miter * (lelt + miter + 1));
  scalar *Z = err + lelt, *M = Z + lelt, *GZ = M + miter * lelt;
  scalar *rhs = GZ + miter * miter, *v = rhs + miter;

  int grammian = 0;

  for (i = 0; i < miter; i++) {
    metric_tic(gsc, PROJECT);
    int ppfi = project(y, lelt, z, wl, d, gsc, 100, 1, 0, buf);
    metric_toc(gsc, PROJECT);
    metric_acc(END, ppfi);

    ortho(y, lelt, nelg, gsc);

    scalar lambda = dot(y, z, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &lambda, 1, bfr);

    for (uint j = 0; j < lelt; j++)
      err[i] = y[i] - lambda * z[i];
    scalar norme = dot(err, err, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &norme, 1, bfr);
    norme = sqrt(norme);

    scalar norm = dot(y, y, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &norm, 1, bfr);
    scalar normi = 1.0 / sqrt(norm);

    for (j = 0; j < lelt; j++)
      z[j] = y[j] * normi;

    ortho(z, lelt, nelg, gsc);

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
            rhs[j] += Z[j * lelt + l] * z[l];
        }
        // Global reduction rhs[j]
        comm_allreduce(gsc, gs_double, gs_add, rhs, i, bfr);

        // Z[k,:] = z[:] - Z[:,1:lelt]*rhs[:]
        for (l = 0; l < lelt; l++)
          Z[i * lelt + l] = z[l];

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
          laplacian(GZ, wl, &Z[j * lelt], buf);
          for (k = 0; k < N; k++) {
            M[k * N + j] = 0.0;
            for (l = 0; l < lelt; l++)
              M[k * N + j] += Z[k * lelt + l] * GZ[l];
          }
        }

        // Global reduction of M
        comm_allreduce(gsc, gs_double, gs_add, M, N * N, buf->ptr);

        // Inverse power iterarion on M
        inv_power_serial(v, N, M, 0);

        for (j = 0; j < lelt; j++)
          z[j] = 0.0;

        for (j = 0; j < N; j++) {
          for (k = 0; k < lelt; k++)
            z[k] += Z[j * lelt + k] * v[j];
        }
        ortho(z, lelt, nelg, gsc);
      } else {
        // Z(k,:) = z;
        for (l = 0; l < lelt; l++)
          Z[i * lelt + l] = z[l];
      }
    }

    if (ppfi == 1)
      break;
  }

  laplacian_free(wl);
  mg_free(d);
  free(err);

  return i == miter ? i : i + 1;
}

static double sign(GenmapScalar a, GenmapScalar b) {
  GenmapScalar m = b >= 0.0 ? 1.0 : -1.0;
  return fabs(a) * m;
}

static int tqli(scalar *eVectors, scalar *eValues, sint n, scalar *diagonal,
                scalar *upper, int id) {
  if (n == 0)
    return 0;

  scalar *d = tcalloc(scalar, 2 * n), *e = d + n;
  sint i;
  for (i = 0; i < n; i++)
    d[i] = diagonal[i];
  for (i = 0; i < n - 1; i++)
    e[i] = upper[i];
  e[n - 1] = 0.0;

  for (i = 0; i < n; i++) {
    for (uint j = 0; j < n; j++)
      eVectors[i * n + j] = 0;
    eVectors[i * n + i] = 1;
  }

  GenmapInt j, k, l, iter, m;
  for (l = 0; l < n; l++) {
    iter = 0;
    do {
      for (m = l; m < n - 1; m++) {
        scalar dd = fabs(d[m]) + fabs(d[m + 1]);
        /* Should use a tolerance for this check */
        if (fabs(e[m]) / dd < GENMAP_TOL)
          break;
      }

      if (m != l) {
        if (iter++ == 30) {
          if (id == 0)
            printf("Too many iterations.\n");
          // vec_copy(*eValues, d);
          for (i = 0; i < n; i++)
            eValues[i] = d[i];
          return 1;
        }

        GenmapScalar g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        GenmapScalar r = sqrt(g * g + 1.0);

        g = d[m] - d[l] + e[l] / (g + sign(r, g));
        GenmapScalar s = 1.0, c = 1.0, p = 0.0;

        for (i = m - 1; i >= l; i--) {
          GenmapScalar f = s * e[i];
          GenmapScalar b = c * e[i];

          if (fabs(f) >= fabs(g)) {
            c = g / f;
            r = sqrt(c * c + 1.0);
            e[i + 1] = f * r;
            s = 1.0 / r;
            c = c * s;
          } else {
            s = f / g;
            r = sqrt(s * s + 1.0);
            e[i + 1] = g * r;
            c = 1.0 / r;
            s = s * c;
          }

          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2.0 * c * b;
          p = s * r;
          d[i + 1] = g + p;
          g = c * r - b;
          /* Find eigenvectors */
          for (k = 0; k < n; k++) {
            f = eVectors[k * n + i + 1];
            eVectors[k * n + i + 1] = s * eVectors[k * n + i] + c * f;
            eVectors[k * n + i] = c * eVectors[k * n + i] - s * f;
          }
          /* Done with eigenvectors */
        }

        if (r < GENMAP_TOL && i >= l)
          continue;

        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while (m != l);
  }

  /* Orthnormalize eigenvectors -- Just normalize? */
  for (i = 0; i < n; i++) {
    for (j = 0; j < i; j++) {
      GenmapScalar tmp = eVectors[i * n + j];
      eVectors[i * n + j] = eVectors[j * n + i];
      eVectors[j * n + i] = tmp;
    }
  }

  for (k = 0; k < n; k++) {
    e[k] = 0;
    for (uint i = 0; i < n; i++)
      e[k] += eVectors[k * n + i] * eVectors[k * n + i];
    if (e[k] > 0.0)
      e[k] = sqrt(fabs(e[k]));
    GenmapScalar scale = 1.0 / e[k];
    for (uint i = 0; i < n; i++)
      eVectors[k * n + i] *= scale;
  }

  // vec_copy(*eValues, d);
  for (i = 0; i < n; i++)
    eValues[i] = d[i];

  free(d);

  return 0;
}

static int lanczos_aux(scalar *diag, scalar *upper, scalar *rr, uint lelt,
                       ulong nelg, int niter, scalar *f, struct laplacian *gl,
                       struct comm *gsc, buffer *bfr) {
  scalar *r = tcalloc(scalar, 3 * lelt), *p = r + lelt, *w = p + lelt;
  // vec_copy(r, f);
  uint i;
  for (i = 0; i < lelt; i++)
    r[i] = f[i];

  // vec_ortho(gsc, r, nelg);
  ortho(r, lelt, nelg, gsc);

  scalar eps = 1e-5, rtz1 = 1, pap = 0, alpha, beta, rtz2, pap_old;

  scalar rtr = dot(r, r, lelt), buf[2];
  comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, buf);
  scalar rnorm = sqrt(rtr), rni = 1.0 / rnorm, rtol = rnorm * eps;

  // vec_scale(rr[0], r, rni);
  for (i = 0; i < lelt; i++)
    rr[0 * lelt + i] = r[i] * rni;

  int iter;
  for (iter = 0; iter < niter; iter++) {
    rtz2 = rtz1, rtz1 = rtr;
    beta = rtz1 / rtz2;
    if (iter == 0)
      beta = 0.0;

    // add2s1(p,r,beta,n)
    for (i = 0; i < lelt; i++)
      p[i] = beta * p[i] + r[i];

    scalar pp = dot(p, p, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &pp, 1, buf);

    // vec_ortho(gsc, p, nelg);
    ortho(p, lelt, nelg, gsc);

    metric_tic(gsc, LAPLACIAN);
    laplacian(w, gl, p, bfr);
    metric_toc(gsc, LAPLACIAN);

    scalar ww = dot(w, w, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &ww, 1, buf);

    pap_old = pap, pap = dot(w, p, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &pap, 1, buf);
#if 0
    if (gsc->id == 0)
      printf("host iter = %d beta = %lf pp = %lf pap = %lf\n", iter, beta, pp,
             pap);
#endif

    alpha = rtz1 / pap;
    // vec_axpby(r, r, 1.0, w, -1.0 * alpha);
    for (i = 0; i < lelt; i++)
      r[i] = r[i] - alpha * w[i];

    rtr = dot(r, r, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, buf);
    rnorm = sqrt(rtr), rni = 1.0 / rnorm;

    // vec_scale(rr[iter + 1], r, rni);
    for (i = 0; i < lelt; i++)
      rr[(iter + 1) * lelt + i] = r[i] * rni;

    if (iter == 0) {
      diag[iter] = pap / rtz1;
    } else {
      diag[iter] = (beta * beta * pap_old + pap) / rtz1;
      upper[iter - 1] = -beta * pap_old / sqrt(rtz2 * rtz1);
    }

    if (rnorm < rtol) {
      iter++;
      break;
    }
  }

  metric_acc(TOL_FNL, rnorm);
  metric_acc(TOL_TGT, rtol);

  free(r);

  return iter;
}

static int lanczos(scalar *fiedler, uint lelt, struct rsb_element *elems,
                   int nv, scalar *initv, struct comm *gsc, int miter,
                   slong nelg, buffer *bfr) {
  metric_tic(gsc, LANCZOS);

  struct laplacian *wl;
#if defined(GENMAP_OCCA)
  wl = laplacian_init(elems, lelt, nv, CSR, gsc, bfr);
  occa_lanczos_init(gsc, wl, miter);
#else
  wl = laplacian_init(elems, lelt, nv, GS, gsc, bfr);
#endif

  if (nelg < miter)
    miter = nelg;

  scalar *alpha = tcalloc(scalar, 2 * miter - 1), *beta = alpha + miter;
  scalar *rr = tcalloc(scalar, (miter + 1) * lelt);
  scalar *eVectors = tcalloc(scalar, miter * miter);
  scalar *eValues = tcalloc(scalar, miter);
  int iter, ipass = 0;
  do {
#if defined(GENMAP_OCCA)
    iter = occa_lanczos_aux(alpha, beta, rr, lelt, nelg, miter, initv, wl, gsc,
                            bfr);
#else
    iter = lanczos_aux(alpha, beta, rr, lelt, nelg, miter, initv, wl, gsc, bfr);
#endif

    // Use TQLI and find the minimum eigenvalue and associated vector
    tqli(eVectors, eValues, iter, alpha, beta, gsc->id);

    GenmapScalar eValMin = fabs(eValues[0]);
    GenmapInt eValMinI = 0;
    uint i;
    for (i = 1; i < iter; i++) {
      if (fabs(eValues[i]) < eValMin) {
        eValMin = fabs(eValues[i]);
        eValMinI = i;
      }
    }

    GenmapInt j;
    for (i = 0; i < lelt; i++) {
      fiedler[i] = 0.0;
      for (j = 0; j < iter; j++)
        fiedler[i] += rr[j * lelt + i] * eVectors[eValMinI * iter + j];
    }

    // vec_ortho(gsc, fiedler, nelg);
    ortho(fiedler, lelt, nelg, gsc);
  } while (iter == miter && ipass++ < miter);

  free(alpha), free(rr), free(eVectors), free(eValues);
#if defined(GENMAP_OCCA)
  occa_lanczos_free();
#endif
  laplacian_free(wl);

  metric_toc(gsc, LANCZOS);

  return ipass;
}

int fiedler(struct rsb_element *elems, uint lelt, int nv, int max_iter,
            int algo, struct comm *gsc, buffer *buf) {
  slong out[2][1], bfr[2][1], in = lelt;
  comm_scan(out, gsc, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0], nelg = out[1][0];

  scalar *initv = tcalloc(scalar, 2 * lelt), *fiedler = initv + lelt;
  uint i;
  for (i = 0; i < lelt; i++) {
    initv[i] = start + i + 1.0;
    // if (start + i < nelg / 2)
    //  initv->data[i] += 1000 * nelg;
  }

  // vec_ortho(gsc, initv, nelg);
  ortho(initv, lelt, nelg, gsc);

  // double rtr = vec_dot(initv, initv);
  scalar rtr = dot(initv, initv, lelt), rni;
  comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, &rni);

  rni = 1.0 / sqrt(rtr);
  // vec_scale(initv, initv, rni);
  for (uint i = 0; i < lelt; i++)
    initv[i] *= rni;

  int iter = 0;
  if (algo == 0)
    iter = lanczos(fiedler, lelt, elems, nv, initv, gsc, max_iter, nelg, buf);
  else if (algo == 1)
    iter = inverse(fiedler, lelt, elems, nv, initv, gsc, max_iter, nelg, buf);
  metric_acc(FIEDLER_NITER, iter);

  GenmapScalar norm = 0;
  for (i = 0; i < lelt; i++)
    norm += fiedler[i] * fiedler[i];
#if 0
  if (gsc->id == 0)
    printf("fiedler norm = %lf\n", norm);
#endif

  GenmapScalar normi;
  comm_allreduce(gsc, gs_double, gs_add, &norm, 1, &normi);
  normi = 1.0 / sqrt(norm);

  // vec_scale(fiedler, fiedler, normi);
  for (i = 0; i < lelt; i++)
    fiedler[i] *= normi;

  for (i = 0; i < lelt; i++)
    elems[i].fiedler = fiedler[i];

  free(initv);
  return 0;
}
