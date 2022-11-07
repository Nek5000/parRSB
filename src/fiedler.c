#include "metrics.h"
#include "multigrid.h"
#include "parrsb-impl.h"
#include "sort.h"

#define MM 500

extern void matrix_inverse(int N, double *A);

inline static scalar dot(scalar *y, scalar *x, uint n) {
  scalar result = 0.0;
  for (uint i = 0; i < n; i++)
    result += x[i] * y[i];

  return result;
}

inline static void ortho(scalar *q, uint n, ulong N, struct comm *c) {
  scalar sum = 0.0;
  for (uint i = 0; i < n; i++)
    sum += q[i];

  scalar wrk;
  comm_allreduce(c, gs_double, gs_add, &sum, 1, &wrk);
  sum /= N;

  for (uint i = 0; i < n; i++)
    q[i] -= sum;
}

inline static void orthol(scalar *q, uint n) {
  scalar sum = 0.0;
  for (uint i = 0; i < n; i++)
    sum += q[i];
  sum /= n;
  for (uint i = 0; i < n; i++)
    q[i] -= sum;
}

struct fiedler {
  scalar fiedler;
  uint proc, seq;
  int part0, part1;
};

int power_serial(double *y, uint N, double *A, int verbose) {
  time_t t;
  srand((unsigned)time(&t));

  int i;
  scalar norm = 0.0;
  for (i = 0; i < N; i++) {
    y[i] = (rand() % 50) / 50.0;
    norm += y[i] * y[i];
  }

  scalar normi = 1.0 / sqrt(norm);
  for (i = 0; i < N; i++)
    y[i] *= normi;

  double *Ay = tcalloc(double, N);
  int j, k, l;
  scalar err = 1.0, lambda;
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
  tfree(Ay);

  return i;
}

int inv_power_serial(double *y, uint N, double *A, int verbose) {
  double *Ainv = tcalloc(double, N *N);
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

  tfree(Ainv);

  return j;
}

static int project(scalar *x, uint n, scalar *b, struct laplacian *L,
                   struct mg *d, struct comm *c, int miter, double tol,
                   int null_space, int verbose, buffer *bfr) {
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
  scalar rtol = rr * tol * tol;

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
    // par_mat_vec(w, p, S, gsh, wrk, bfr);
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

    mg_vcycle(z, r, d, c, bfr);

    rzt = rz1;
    if (null_space)
      ortho(z, n, ng, c);
    rz1 = dot(r, z, n);
    comm_allreduce(c, gs_double, gs_add, &rz1, 1, buf);

    for (j = 0; j < n; j++)
      dz[j] = z[j] - z0[j];
    rz2 = dot(r, dz, n);
    comm_allreduce(c, gs_double, gs_add, &rz2, 1, buf);

    if (c->id == 0 && verbose > 0) {
      printf("rr = %lf rtol = %lf rz0 = %lf rz1 = %lf rz2 = %lf\n", rr, rtol,
             rzt, rz1, rz2);
      fflush(stdout);
    }

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

  tfree(z), tfree(P);

  return i == miter ? i : i + 1;
}

static struct par_mat *setup_csr_from_con(const uint nelt, const ulong *eid,
                                          const slong *vtx, int nv, int sep,
                                          struct comm *c, struct crystal *cr,
                                          buffer *bfr) {
  struct array nbrs, eij;
  array_init(struct nbr, &nbrs, 100);
  array_init(struct mij, &eij, 100);

  find_nbrs(&nbrs, eid, vtx, nelt, nv, cr, bfr);
  compress_nbrs(&eij, &nbrs, bfr);
  array_free(&nbrs);

  struct par_mat *M = tcalloc(struct par_mat, 1);
  par_csr_setup(M, &eij, sep, bfr);
  array_free(&eij);

  return M;
}

// Input z should be orthogonal to 1-vector, have unit norm.
// inverse iteration should not change z.
static int inverse(scalar *y, struct array *elements, int nv, scalar *z,
                   struct comm *gsc, int miter, int mpass, double tol,
                   int factor, int sagg, int grammian, slong nelg,
                   buffer *buf) {
  metric_tic(gsc, RSB_INVERSE_SETUP);
  uint lelt = elements->n;
  struct rsb_element *elems = (struct rsb_element *)elements->ptr;
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

  // Setup LAMG preconditioner
  struct crystal cr;
  crystal_init(&cr, gsc);
  struct par_mat *L = setup_csr_from_con(lelt, eid, vtx, nv, 1, gsc, &cr, buf);
  for (uint i = 0; i < L->rn; i++) {
    for (uint j = L->adj_off[i], je = L->adj_off[i + 1]; j < je; j++)
      L->adj_val[j] /= L->diag_val[i];
    L->diag_val[i] = 1.0;
  }
  struct mg *d = mg_setup(L, factor, sagg, &cr, buf);
  crystal_free(&cr);
  metric_toc(gsc, RSB_INVERSE_SETUP);

  scalar *err = tcalloc(scalar, 2 * lelt + miter * (lelt + miter + 1));
  scalar *Z = err + lelt, *M = Z + lelt, *GZ = M + miter * lelt;
  scalar *rhs = GZ + miter * miter, *v = rhs + miter;

  metric_tic(gsc, RSB_INVERSE);
  for (i = 0; i < mpass; i++) {
    int ppfi = project(y, lelt, z, wl, d, gsc, miter, tol, 1, 0, buf);

    ortho(y, lelt, nelg, gsc);

    scalar lambda = dot(y, z, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &lambda, 1, bfr);

    for (uint j = 0; j < lelt; j++)
      err[j] = y[j] - lambda * z[j];
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
  metric_toc(gsc, RSB_INVERSE);

  if (L) {
    par_mat_free(L);
    tfree(L);
  }
  laplacian_free(wl), mg_free(d), tfree(err);

  return i;
}

static double sign(scalar a, scalar b) {
  scalar m = b >= 0.0 ? 1.0 : -1.0;
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

  int j, k, l, iter, m;
  for (l = 0; l < n; l++) {
    iter = 0;
    do {
      for (m = l; m < n - 1; m++) {
        scalar dd = fabs(d[m]) + fabs(d[m + 1]);
        /* Should use a tolerance for this check */
        if (fabs(e[m]) / dd < SCALAR_TOL)
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

        scalar g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        scalar r = sqrt(g * g + 1.0);

        g = d[m] - d[l] + e[l] / (g + sign(r, g));
        scalar s = 1.0, c = 1.0, p = 0.0;

        for (i = m - 1; i >= l; i--) {
          scalar f = s * e[i];
          scalar b = c * e[i];

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

        if (r < SCALAR_TOL && i >= l)
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
      scalar tmp = eVectors[i * n + j];
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
    scalar scale = 1.0 / e[k];
    for (uint i = 0; i < n; i++)
      eVectors[k * n + i] *= scale;
  }

  // vec_copy(*eValues, d);
  for (i = 0; i < n; i++)
    eValues[i] = d[i];

  tfree(d);

  return 0;
}

static int lanczos_aux(scalar *diag, scalar *upper, scalar *rr, uint n,
                       ulong ng, int niter, double tol, scalar *f,
                       struct laplacian *gl, struct comm *gsc, buffer *bfr) {
  scalar *r = tcalloc(scalar, n);
  scalar *p = tcalloc(scalar, n);
  scalar *w = tcalloc(scalar, n);

  for (uint i = 0; i < n; i++)
    r[i] = f[i];

  ortho(r, n, ng, gsc);

  scalar rtz1 = 1, rtz2, pap1 = 0, pap2;
  scalar rtr = dot(r, r, n), wrk[2];
  comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, wrk);
  scalar rnorm = sqrt(rtr), rtol = rnorm * tol;

  // vec_scale(rr[0], r, rni);
  scalar rni = 1.0 / rnorm;
  for (uint i = 0; i < n; i++)
    rr[0 * n + i] = r[i] * rni;

  unsigned iter = 0;
  while (iter < niter) {
    rtz2 = rtz1, rtz1 = rtr;
    scalar beta = rtz1 / rtz2;
    if (iter == 0)
      beta = 0.0;

    // add2s1(p, r, beta, n)
    for (uint i = 0; i < n; i++)
      p[i] = beta * p[i] + r[i];

    // scalar pp = dot(p, p, n);
    // comm_allreduce(gsc, gs_double, gs_add, &pp, 1, wrk);

    // vec_ortho(gsc, p, ng);
    ortho(p, n, ng, gsc);

    laplacian(w, gl, p, bfr);

    // scalar ww = dot(w, w, n);
    // comm_allreduce(gsc, gs_double, gs_add, &ww, 1, wrk);

    pap2 = pap1, pap1 = dot(w, p, n);
    comm_allreduce(gsc, gs_double, gs_add, &pap1, 1, wrk);

    // vec_axpby(r, r, 1.0, w, -1.0 * alpha);
    scalar alpha = rtz1 / pap1;
    for (uint i = 0; i < n; i++)
      r[i] = r[i] - alpha * w[i];

    rtr = dot(r, r, n);
    comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, wrk);
    rnorm = sqrt(rtr);

    // vec_scale(rr[iter + 1], r, rni);
    scalar rni = 1.0 / rnorm;
    for (uint i = 0; i < n; i++)
      rr[(iter + 1) * n + i] = r[i] * rni;

    if (iter == 0) {
      diag[iter] = pap1 / rtz1;
    } else {
      diag[iter] = (beta * beta * pap2 + pap1) / rtz1;
      upper[iter - 1] = -beta * pap2 / sqrt(rtz2 * rtz1);
    }

    iter++;
    if (rnorm < rtol)
      break;
  }

  metric_acc(TOL_FNL, rnorm);
  metric_acc(TOL_TGT, rtol);

  tfree(r), tfree(p), tfree(w);
  return iter;
}

static int lanczos(scalar *fiedler, struct array *arr, int nv, scalar *initv,
                   struct comm *c, int miter, int mpass, double tol, slong ng,
                   buffer *bfr) {
  metric_tic(c, RSB_LANCZOS_SETUP);
  uint n = arr->n;
  struct rsb_element *pa = (struct rsb_element *)arr->ptr;
  struct laplacian *wl = laplacian_init(pa, n, nv, GS, c, bfr);
  metric_toc(c, RSB_LANCZOS_SETUP);

  if (miter > ng)
    miter = ng;

  scalar *alpha = tcalloc(scalar, miter);
  scalar *beta = tcalloc(scalar, miter);
  scalar *rr = tcalloc(scalar, (miter + 1) * n);
  scalar *evecs = tcalloc(scalar, miter * miter);
  scalar *evals = tcalloc(scalar, miter);

  if (mpass < 1)
    mpass = 1;

  unsigned iter = miter, pass = 0;
  while (pass < mpass) {
    double t = comm_time();
    iter = lanczos_aux(alpha, beta, rr, n, ng, miter, tol, initv, wl, c, bfr);
    metric_acc(RSB_LANCZOS, comm_time() - t);

    // Use TQLI and find the minimum eigenvalue and associated vector
    t = comm_time();
    tqli(evecs, evals, iter, alpha, beta, c->id);

    scalar ev_min = fabs(evals[0]);
    uint ev_min_idx = 0;
    for (uint i = 1; i < iter; i++) {
      if (fabs(evals[i]) < ev_min) {
        ev_min = fabs(evals[i]);
        ev_min_idx = i;
      }
    }

    for (uint i = 0; i < n; i++) {
      fiedler[i] = 0.0;
      for (uint j = 0; j < iter; j++)
        fiedler[i] += rr[j * n + i] * evecs[ev_min_idx * iter + j];
    }
    ortho(fiedler, n, ng, c);

    pass++;
    if (iter < miter)
      break;

    for (uint i = 0; i < n; i++)
      initv[i] = fiedler[i];

    metric_acc(RSB_LANCZOS_TQLI, comm_time() - t);
  }

  tfree(alpha), tfree(beta), tfree(rr), tfree(evecs), tfree(evals);
  laplacian_free(wl);

  return (pass - 1) * miter + iter;
}

void fiedler(struct array *arr, int nv, parrsb_options *opts, struct comm *c,
             buffer *bfr, int verbose) {
  metric_tic(c, RSB_FIEDLER_SETUP);
  slong out[2][1], wrk[2][1], in = arr->n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0], ng = out[1][0];

  uint n = arr->n;
  scalar *initv = tcalloc(scalar, n);
  for (uint i = 0; i < n; i++) {
    initv[i] = start + i + 1.0;
    if (start + i < ng / 2)
      initv[i] += 1000 * ng;
  }

  ortho(initv, n, ng, c);
  scalar rtr = dot(initv, initv, n);
  comm_allreduce(c, gs_double, gs_add, &rtr, 1, wrk);

  scalar rni = 1.0 / sqrt(rtr);
  for (uint i = 0; i < n; i++)
    initv[i] *= rni;
  metric_toc(c, RSB_FIEDLER_SETUP);

  metric_tic(c, RSB_FIEDLER_CALC);
  unsigned iter = 0;
  scalar *f = tcalloc(scalar, n);
  switch (opts->rsb_algo) {
  case 0:
    iter = lanczos(f, arr, nv, initv, c, opts->rsb_max_iter,
                   opts->rsb_max_passes, opts->rsb_tol, ng, bfr);
    break;
  case 1:
    iter = inverse(f, arr, nv, initv, c, opts->rsb_max_iter,
                   opts->rsb_max_passes, opts->rsb_tol, opts->rsb_mg_factor,
                   opts->rsb_mg_sagg, opts->rsb_mg_grammian, ng, bfr);
    break;
  default:
    break;
  }
  metric_toc(c, RSB_FIEDLER_CALC);
  metric_acc(RSB_FIEDLER_CALC_NITER, iter);

  scalar norm = dot(f, f, n);
  comm_allreduce(c, gs_double, gs_add, &norm, 1, wrk);

  rni = 1.0 / sqrt(norm);
  for (uint i = 0; i < n; i++)
    f[i] *= rni;

  struct rsb_element *pa = (struct rsb_element *)arr->ptr;
  for (uint i = 0; i < n; i++)
    pa[i].fiedler = f[i];

  tfree(initv), tfree(f);
}

static struct mat *laplacian_local(const struct rsb_element *pe, uint ne,
                                   unsigned nv, unsigned sep, buffer *bfr) {
  struct nbr_t {
    ulong r, c;
  };

  struct array vertices;
  array_init(struct nbr_t, &vertices, ne * nv + 1);

  struct nbr_t v = {.r = 0, .c = 0};
  for (uint i = 0; i < ne; i++) {
    v.r = pe[i].globalId;
    for (unsigned j = 0; j < nv; j++) {
      v.c = pe[i].vertices[j];
      array_cat(struct nbr_t, &vertices, &v, 1);
    }
  }

  struct array nbrs;
  array_init(struct nbr_t, &nbrs, vertices.n * 10 + 1);

  sarray_sort(struct nbr_t, vertices.ptr, vertices.n, c, 1, bfr);
  const struct nbr_t *pv = (struct nbr_t *)vertices.ptr;

  sint i = 0;
  while (i < vertices.n) {
    sint j = i + 1;
    while (j < vertices.n && pv[j].c == pv[i].c)
      j++;
    for (unsigned s = i; s < j; s++) {
      for (unsigned t = i; t < j; t++) {
        v.r = pv[s].r, v.c = pv[t].r;
        array_cat(struct nbr_t, &nbrs, &v, 1);
      }
    }
  }
  array_free(&vertices);

  struct array mijs;
  array_init(struct mij, &mijs, nbrs.n * 10 + 1);

  sarray_sort_2(struct nbr_t, nbrs.ptr, nbrs.n, r, 1, c, 1, bfr);
  const struct nbr_t *pn = (struct nbr_t *)nbrs.ptr;

  struct mij m = {.r = 0, .c = 0, .v = 0};
  i = 0;
  while (i < nbrs.n) {
    sint j = i + 1;
    while (j < nbrs.n && pn[i].c == pn[j].c && pn[i].r == pn[j].r)
      j++;
    m.r = pn[i].r, pn[i].c, m.v = i - j; // -(j - i)
    array_cat(struct mij, &mijs, &m, 1);
  }
  array_free(&nbrs);

  // Make sure that we have zero row sum
  struct mij *pm = (struct mij *)mijs.ptr;
  i = 0;
  while (i < mijs.n) {
    sint j = i, d = -1, sum = 0;
    // Go through the row
    while (j < mijs.n && pm[j].r == pm[i].r) {
      if (pm[j].r == pm[j].c)
        d = j;
      else
        sum += pm[j].v;
      j++;
    }
    assert(d >= 0);
    pm[d].v = -sum;
    i = j;
  }

  struct mat *A = tcalloc(struct mat, 1);
  mat_setup(A, &mijs, 0, bfr);

  array_free(&mijs);

  return A;
}

static int lanczos_aux_local(scalar *diag, scalar *upper, scalar *rr, uint n,
                             int niter, double tol, scalar *f, struct mat *L,
                             buffer *bfr) {
  scalar *r = tcalloc(scalar, n);
  scalar *p = tcalloc(scalar, n);
  scalar *w = tcalloc(scalar, n);

  for (uint i = 0; i < n; i++)
    r[i] = f[i];

  orthol(r, n);

  scalar rtz1 = 1, rtz2, pap1 = 0, pap2;
  scalar rtr = dot(r, r, n), wrk[2];
  scalar rnorm = sqrt(rtr), rtol = rnorm * tol;

  // vec_scale(rr[0], r, rni);
  scalar rni = 1.0 / rnorm;
  for (uint i = 0; i < n; i++)
    rr[0 * n + i] = r[i] * rni;

  unsigned iter = 0;
  while (iter < niter) {
    rtz2 = rtz1, rtz1 = rtr;
    scalar beta = rtz1 / rtz2;
    if (iter == 0)
      beta = 0.0;

    // add2s1(p, r, beta, n)
    for (uint i = 0; i < n; i++)
      p[i] = beta * p[i] + r[i];

    // scalar pp = dot(p, p, n);

    // vec_ortho(gsc, p, ng);
    orthol(p, n);

    mat_vec(w, L, p);

    // scalar ww = dot(w, w, n);

    pap2 = pap1, pap1 = dot(w, p, n);

    // vec_axpby(r, r, 1.0, w, -1.0 * alpha);
    scalar alpha = rtz1 / pap1;
    for (uint i = 0; i < n; i++)
      r[i] = r[i] - alpha * w[i];

    rtr = dot(r, r, n), rnorm = sqrt(rtr);

    // vec_scale(rr[iter + 1], r, rni);
    scalar rni = 1.0 / rnorm;
    for (uint i = 0; i < n; i++)
      rr[(iter + 1) * n + i] = r[i] * rni;

    if (iter == 0) {
      diag[iter] = pap1 / rtz1;
    } else {
      diag[iter] = (beta * beta * pap2 + pap1) / rtz1;
      upper[iter - 1] = -beta * pap2 / sqrt(rtz2 * rtz1);
    }

    iter++;
    if (rnorm < rtol)
      break;
  }

  tfree(r), tfree(p), tfree(w);
  return iter;
}

static int lanczos_local(scalar *fiedler, struct array *arr, uint sidx,
                         uint eidx, int nv, scalar *initv, int miter, int mpass,
                         double tol, buffer *bfr) {
  if (eidx <= sidx + 1)
    return 0;

  uint n = eidx - sidx;
  struct rsb_element *pa = (struct rsb_element *)arr->ptr;
  struct mat *A = laplacian_local(&pa[sidx], n, nv, 0, bfr);

  if (miter > n)
    miter = n;

  scalar *alpha = tcalloc(scalar, miter);
  scalar *beta = tcalloc(scalar, miter);
  scalar *rr = tcalloc(scalar, (miter + 1) * n);
  scalar *evecs = tcalloc(scalar, miter * miter);
  scalar *evals = tcalloc(scalar, miter);

  if (mpass < 1)
    mpass = 1;

  unsigned iter = miter, pass = 0;
  while (pass < mpass) {
    double t = comm_time();
    iter = lanczos_aux_local(alpha, beta, rr, n, miter, tol, initv, A, bfr);

    // Use TQLI and find the minimum eigenvalue and associated vector
    tqli(evecs, evals, iter, alpha, beta, 0);

    scalar ev_min = fabs(evals[0]);
    uint ev_min_idx = 0;
    for (uint i = 1; i < iter; i++) {
      if (fabs(evals[i]) < ev_min) {
        ev_min = fabs(evals[i]);
        ev_min_idx = i;
      }
    }

    for (uint i = 0; i < n; i++) {
      fiedler[i] = 0.0;
      for (uint j = 0; j < iter; j++)
        fiedler[i] += rr[j * n + i] * evecs[ev_min_idx * iter + j];
    }
    orthol(fiedler, n);

    pass++;
    if (iter < miter)
      break;

    for (uint i = 0; i < n; i++)
      initv[i] = fiedler[i];
  }

  tfree(alpha), tfree(beta), tfree(rr), tfree(evecs), tfree(evals);
  mat_free(A);

  return (pass - 1) * miter + iter;
}

void fiedler_local(struct array *arr, uint sidx, uint eidx, int nv,
                   parrsb_options *opts, buffer *bfr, int verbose) {
  if (eidx <= sidx + 1)
    return;

  uint n = eidx - sidx;
  scalar *initv = tcalloc(scalar, n);
  for (uint i = 0; i < n / 2; i++)
    initv[i] = i + 1.0;
  for (uint i = n / 2; i < n; i++)
    initv[i] = i + 1.0 + 1000 * n;

  orthol(initv, n);
  scalar rtr = dot(initv, initv, n);
  rtr = 1.0 / sqrt(rtr);
  for (uint i = 0; i < n; i++)
    initv[i] *= rtr;

  unsigned iter = 0;
  scalar *f = tcalloc(scalar, n);
  switch (opts->rsb_algo) {
  case 0:
    iter = lanczos_local(f, arr, sidx, eidx, nv, initv, opts->rsb_max_iter,
                         opts->rsb_max_passes, opts->rsb_tol, bfr);
    break;
  default:
    break;
  }

  scalar norm = 1.0 / sqrt(dot(f, f, n));
  for (uint i = 0; i < n; i++)
    f[i] *= norm;

  struct rsb_element *pa = (struct rsb_element *)arr->ptr;
  for (uint i = 0; i < n; i++)
    pa[sidx + i].fiedler = f[i];

  tfree(initv), tfree(f);

  uint midx = (sidx + eidx) / 2;
  fiedler_local(arr, sidx, midx, nv, opts, bfr, verbose);
  fiedler_local(arr, midx, eidx, nv, opts, bfr, verbose);
}
