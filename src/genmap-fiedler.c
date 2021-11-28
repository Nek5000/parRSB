#include <math.h>
#include <time.h>

#include <genmap-impl.h>
#include <genmap-multigrid.h>
#include <sort.h>

#define MM 500

static int vec_create(genmap_vector *x, GenmapInt size) {
  /* Asserts:
       - size > 0
  */
  assert(size > 0);

  GenmapMalloc(1, x);
  if (*x == NULL) {
    return 1;
  }

  (*x)->size = size;
  (*x)->data = NULL;

  GenmapMalloc((size_t)size, &(*x)->data);
  if ((*x)->data == NULL) {
    return 1;
  }

  return 0;
}

static int vec_destroy(genmap_vector x) {
  if (x->data) {
    GenmapFree(x->data);
  }

  if (x) {
    GenmapFree(x);
  }

  return 0;
}

static int vec_create_zeros(genmap_vector *x, GenmapInt size) {
  vec_create(x, size);

  GenmapInt i;
  for (i = 0; i < size; i++) {
    (*x)->data[i] = 0.;
  }

  return 0;
}

static int vec_copy(genmap_vector y, genmap_vector x) {
  /* Asserts:
       - size y = size x
  */
  assert(y->size >= x->size);

  GenmapInt n = x->size;
  GenmapInt i;
  for (i = 0; i < n; i++) {
    y->data[i] = x->data[i];
  }

  return 0;
}

static int vec_scale(genmap_vector y, genmap_vector x, GenmapScalar alpha) {
  /* asserts:
       - size x = size y
  */
  assert(x->size == y->size);

  GenmapInt n = x->size;
  GenmapInt i;
  for (i = 0; i < n; i++) {
    y->data[i] = alpha * x->data[i];
  }

  return 0;
}

static GenmapScalar vec_dot(genmap_vector y, genmap_vector x) {
  /* asserts:
       - size x = size y
  */
  assert(x->size == y->size);

  GenmapScalar result = 0.0;
  GenmapInt i;
  for (i = 0; i < x->size; i++) {
    result += x->data[i] * y->data[i];
  }

  return result;
}

static int vec_axpby(genmap_vector z, genmap_vector x, GenmapScalar alpha,
                     genmap_vector y, GenmapScalar beta) {
  assert(z->size == x->size);
  assert(z->size == y->size);

  GenmapInt n = z->size;
  GenmapInt i;
  for (i = 0; i < n; i++) {
    z->data[i] = alpha * x->data[i] + beta * y->data[i];
  }

  return 0;
}

/* Orthogonalize by 1-vector (vector of all 1's) */
static int vec_ortho(struct comm *c, genmap_vector q1, GenmapULong n) {
  GenmapInt i;
  GenmapScalar sum = 0.0;
  for (i = 0; i < q1->size; i++)
    sum += q1->data[i];

  GenmapScalar buf;
  comm_allreduce(c, gs_double, gs_add, &sum, 1, &buf);
  sum /= n;

  for (i = 0; i < q1->size; i++)
    q1->data[i] -= sum;

  return 0;
}

static int vec_dssum(GenmapScalar *v, struct laplacian *gl, buffer *buf) {
  uint lelt = gl->nel;
  int nv = gl->nv;

  GenmapInt i, j;
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      gl->u[nv * i + j] = v[i];

  gs(gl->u, gs_double, gs_add, 0, gl->gsh, buf);

  for (i = 0; i < lelt; i++) {
    v[i] = 0.0;
    for (j = 0; j < nv; j++)
      v[i] += gl->u[nv * i + j];
  }

  return 0;
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

  slong out[2][1], bfr[2][1];
  slong in = fdlr->n;
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

int project(genmap_vector x, struct laplacian *gl, struct mg_data *d,
            genmap_vector ri, int max_iter, struct comm *gsc, buffer *buf,
            int gid) {
  assert(x->size == ri->size);
  uint lelt = x->size;

  genmap_vector z0, z, dz, w, p, r;
  vec_create(&z, lelt);
  vec_create(&w, lelt);
  vec_create(&r, lelt);
  vec_create(&p, lelt);
  vec_create(&z0, lelt);
  vec_create(&dz, lelt);

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

  GenmapScalar rr = vec_dot(r, r);
  comm_allreduce(gsc, gs_double, gs_add, &rr, 1, bfr);
  double tol = 1e-3;
  double res_tol = rr * tol * tol;

  vec_copy(z, r);
  vec_ortho(gsc, z, nelg);

  GenmapScalar rz1 = vec_dot(r, z);
  comm_allreduce(gsc, gs_double, gs_add, &rz1, 1, bfr);

  vec_copy(p, z);

  GenmapScalar alpha, beta, rz0, rz2, scale;

  uint j, k;
  i = 0;
  while (i < max_iter) {
    metric_tic(gsc, LAPLACIAN);
    laplacian(w->data, gl, p->data, buf);
    metric_toc(gsc, LAPLACIAN);

    GenmapScalar den = vec_dot(p, w);
    comm_allreduce(gsc, gs_double, gs_add, &den, 1, bfr);
    alpha = rz1 / den;

    scale = 1.0 / sqrt(den);
    for (j = 0; j < lelt; j++) {
      W[i * lelt + j] = scale * w->data[j];
      P[i * lelt + j] = scale * p->data[j];
    }

    vec_axpby(x, x, 1.0, p, alpha);
    vec_axpby(r, r, 1.0, w, -alpha);

    rr = vec_dot(r, r);
    comm_allreduce(gsc, gs_double, gs_add, &rr, 1, bfr);

    if (rr < res_tol || sqrt(rr) < tol)
      break;

    GenmapScalar norm0 = vec_dot(z, z);
    comm_allreduce(gsc, gs_double, gs_add, &norm0, 1, bfr);

    vec_copy(z0, z);
#if 1
    mg_vcycle(z->data, r->data, d);
#else
    vec_copy(z, r);
#endif

    GenmapScalar norm1 = vec_dot(z, z);
    comm_allreduce(gsc, gs_double, gs_add, &norm1, 1, bfr);

    rz0 = rz1;
    vec_ortho(gsc, z, nelg);
    rz1 = vec_dot(r, z);
    comm_allreduce(gsc, gs_double, gs_add, &rz1, 1, bfr);

    vec_axpby(dz, z, 1.0, z0, -1.0);
    rz2 = vec_dot(r, dz);
    comm_allreduce(gsc, gs_double, gs_add, &rz2, 1, bfr);

    beta = rz2 / rz0;
    vec_axpby(p, z, 1.0, p, beta);

    GenmapScalar normw = vec_dot(w, w);
    comm_allreduce(gsc, gs_double, gs_add, &normw, 1, bfr);

    i++;

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
  }

  vec_destroy(z);
  vec_destroy(w);
  vec_destroy(p);
  vec_destroy(r);
  vec_destroy(z0);
  vec_destroy(dz);

  GenmapFree(P);
  GenmapFree(W);

  return i == max_iter ? i : i + 1;
}

// Input z should be orthogonal to 1-vector, have unit norm.
// inverse iteration should not change z.
static int inverse(genmap_vector y, struct laplacian *gl, struct mg_data *d,
                   genmap_vector z, struct comm *gsc, int max_iter,
                   int grammian, struct array *fdlr, slong nelg, buffer *bff,
                   int gid) {
  assert(z->size == y->size);
  assert(z->size == gl->nel);

  uint lelt = gl->nel;

  // Grammian
  GenmapScalar *Z, *GZ, *M, *rhs, *v, *buf;
  GenmapMalloc(max_iter * lelt, &Z);
  GenmapMalloc(lelt, &GZ);
  GenmapMalloc(max_iter * max_iter, &M);
  GenmapMalloc(max_iter, &rhs);
  GenmapMalloc(max_iter, &v);
  GenmapMalloc(max_iter * max_iter, &buf);

  genmap_vector err;
  vec_create(&err, lelt);

  slong bfr[2];
  uint i, j, k, l;
  for (i = 0; i < max_iter; i++) {
    metric_tic(gsc, PROJECT);
    int ppfi = project(y, gl, d, z, 100, gsc, bff, gid);
    metric_toc(gsc, PROJECT);

    vec_ortho(gsc, y, nelg);

    GenmapScalar lambda = vec_dot(y, z);
    comm_allreduce(gsc, gs_double, gs_add, &lambda, 1, bfr);

    vec_axpby(err, y, 1.0, z, -lambda);
    GenmapScalar norme = vec_dot(err, err);
    comm_allreduce(gsc, gs_double, gs_add, &norme, 1, bfr);
    norme = sqrt(norme);

    metric_acc(END, ppfi);

    GenmapScalar norm = vec_dot(y, y);
    comm_allreduce(gsc, gs_double, gs_add, &norm, 1, bfr);
    GenmapScalar normi = 1.0 / sqrt(norm);

    for (j = 0; j < lelt; j++)
      z->data[j] = y->data[j] * normi;

    vec_ortho(gsc, z, nelg);

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
          laplacian(GZ, gl, &Z[j * lelt], bff);
          for (k = 0; k < N; k++) {
            M[k * N + j] = 0.0;
            for (l = 0; l < lelt; l++)
              M[k * N + j] += Z[k * lelt + l] * GZ[l];
          }
        }

        // Global reduction of M
        comm_allreduce(gsc, gs_double, gs_add, M, N * N, buf);

        // Inverse power iterarion on M
        inv_power_serial(v, N, M, 0);

        for (j = 0; j < lelt; j++)
          z->data[j] = 0.0;

        for (j = 0; j < N; j++) {
          for (k = 0; k < lelt; k++)
            z->data[k] += Z[j * lelt + k] * v[j];
        }
        vec_ortho(gsc, z, nelg);
      } else {
        // Z(k,:) = z;
        for (l = 0; l < lelt; l++)
          Z[i * lelt + l] = z->data[l];
      }
    }

    if (ppfi == 1)
      break;
  }

  GenmapFree(Z);
  GenmapFree(GZ);
  GenmapFree(M);
  GenmapFree(rhs);
  GenmapFree(v);
  GenmapFree(buf);

  vec_destroy(err);

  return i == max_iter ? i : i + 1;
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
  vec_create(&d, n);
  vec_create(&e, n);

  vec_copy(d, diagonal);
  vec_copy(e, upper);
  e->data[n - 1] = 0.0;

  /* Create the vector to store eigenvalues */
  vec_create(eValues, n);
  /* Init to identity */
  GenmapMalloc(n, eVectors);
  GenmapInt i;
  for (i = 0; i < n; i++) {
    vec_create_zeros(&(*eVectors)[i], n);
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
          vec_copy(*eValues, d);
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
    e->data[k] = vec_dot((*eVectors)[k], (*eVectors)[k]);
    if (e->data[k] > 0.0)
      e->data[k] = sqrt(fabs(e->data[k]));
    GenmapScalar scale = 1.0 / e->data[k];
    vec_scale((*eVectors)[k], (*eVectors)[k], scale);
  }

  vec_copy(*eValues, d);

  vec_destroy(d);
  vec_destroy(e);

  return 0;
}

static int lanczos_aux(genmap_vector diag, genmap_vector upper,
                       genmap_vector *rr, uint lelt, ulong nelg, int niter,
                       genmap_vector f, struct laplacian *gl, struct comm *gsc,
                       buffer *bfr) {
  assert(f->size == lelt);

  genmap_vector r, p, w;
  vec_create_zeros(&p, lelt);
  vec_create(&w, lelt);
  vec_create(&r, lelt);

  vec_copy(r, f);

  vec_ortho(gsc, r, nelg);

  GenmapScalar buf[2];
  GenmapScalar rtr = vec_dot(r, r);
  comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, buf);
  GenmapScalar rnorm = sqrt(rtr);
  GenmapScalar rni = 1.0 / rnorm;
  vec_scale(rr[0], r, rni);

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
    uint i;
    for (i = 0; i < lelt; i++)
      p->data[i] = beta * p->data[i] + r->data[i];

    GenmapScalar pp = vec_dot(p, p);
    comm_allreduce(gsc, gs_double, gs_add, &pp, 1, buf);

    vec_ortho(gsc, p, nelg);

    laplacian(w->data, gl, p->data, bfr);

    GenmapScalar ww = vec_dot(w, w);
    comm_allreduce(gsc, gs_double, gs_add, &ww, 1, buf);

    pap_old = pap;
    pap = vec_dot(w, p);
    comm_allreduce(gsc, gs_double, gs_add, &pap, 1, buf);
#if 0
    if (gsc->id == 0)
      printf("iter = %d beta = %lf pp = %lf pap = %lf ww = %lf\n", iter,
             beta, pp, pap, ww);
#endif

    alpha = rtz1 / pap;
    vec_axpby(r, r, 1.0, w, -1.0 * alpha);

    rtr = vec_dot(r, r);
    comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, buf);
    rnorm = sqrt(rtr);
    rni = 1.0 / rnorm;

    vec_scale(rr[iter + 1], r, rni);

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

  vec_destroy(p);
  vec_destroy(w);
  vec_destroy(r);

  return iter;
}

static int lanczos(genmap_vector fiedler, struct laplacian *gl,
                   genmap_vector initv, struct comm *gsc, int max_iter,
                   struct array *fdlr, slong nelg, buffer *bfr) {
  metric_tic(gsc, LANCZOS);

  genmap_vector alpha, beta;
  vec_create(&alpha, max_iter);
  vec_create(&beta, max_iter - 1);

  if (nelg < max_iter) {
    max_iter = nelg;
    alpha->size = max_iter;
    beta->size = max_iter - 1;
  }

  genmap_vector *rr;
  GenmapMalloc(max_iter + 1, &rr);

  uint lelt = initv->size;
  uint i;
  for (i = 0; i < max_iter + 1; ++i)
    vec_create(&rr[i], lelt);

  int iter, ipass = 0;
  do {
    iter =
        lanczos_aux(alpha, beta, rr, lelt, nelg, max_iter, initv, gl, gsc, bfr);

    genmap_vector evTriDiag;
    vec_create(&evTriDiag, iter);

    /* Use TQLI and find the minimum eigenvalue and associated vector */
    genmap_vector *eVectors, eValues;
    tqli(&eVectors, &eValues, alpha, beta, gsc->id);

    GenmapScalar eValMin = fabs(eValues->data[0]);
    GenmapInt eValMinI = 0;
    uint i;
    for (i = 1; i < iter; i++) {
      if (fabs(eValues->data[i]) < eValMin) {
        eValMin = fabs(eValues->data[i]);
        eValMinI = i;
      }
    }
    vec_copy(evTriDiag, eVectors[eValMinI]);

    GenmapInt j;
    for (i = 0; i < lelt; i++) {
      fiedler->data[i] = 0.0;
      for (j = 0; j < iter; j++)
        fiedler->data[i] += rr[j]->data[i] * evTriDiag->data[j];
    }

    vec_ortho(gsc, fiedler, nelg);

    vec_destroy(evTriDiag);
    vec_destroy(eValues);
    for (i = 0; i < iter; i++)
      vec_destroy(eVectors[i]);
    GenmapFree(eVectors);
  } while (iter == max_iter && ipass++ < max_iter);

  metric_toc(gsc, LANCZOS);

  for (i = 0; i < max_iter + 1; ++i)
    vec_destroy(rr[i]);
  vec_destroy(alpha);
  vec_destroy(beta);

  return ipass;
}

int fiedler(struct rsb_element *elems, uint lelt, int nv, int max_iter,
            int algo, struct comm *gsc, buffer *buf, int gid) {
  slong out[2][1], bfr[2][1];
  slong in = lelt;
  comm_scan(out, gsc, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0];
  slong nelg = out[1][0];

  struct laplacian l;
  struct mg_data d;
  if (algo == 1) {
    laplacian_init(&l, elems, lelt, nv, CSR | UNWEIGHTED, gsc, buf);
    mg_setup(&d, 4, gsc, l.M);
  }

  struct laplacian wl;
#if 0
  laplacian_init(&wl, elems, lelt, nv, GS | WEIGHTED, gsc, buf);
#else
  // struct laplacian gsl;
  // laplacian_init(&gsl, elems, lelt, nv, GS | WEIGHTED, gsc, buf);
  // laplacian_free(&gsl);
  laplacian_init(&wl, elems, lelt, nv, CSR | WEIGHTED, gsc, buf);
  // csr_mat_print(wl.M, gsc);
#endif

#if defined(GENMAP_OCCA)
  occa_lanczos_init(gsc, &wl, max_iter);
#endif

  genmap_vector initv, fiedler;
  vec_create(&initv, lelt);
  vec_create(&fiedler, lelt);

  uint i;
  for (i = 0; i < lelt; i++) {
    initv->data[i] = start + i + 1.0;
    // if (start + i < nelg / 2)
    //  initv->data[i] += 1000 * nelg;
  }

  vec_ortho(gsc, initv, nelg);

  double rtr = vec_dot(initv, initv);
  double rni;
  comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, &rni);

  rni = 1.0 / sqrt(rtr);
  vec_scale(initv, initv, rni);

  struct array fdlr;
  array_init(struct fiedler, &fdlr, lelt);

  fdlr.n = lelt;
  struct fiedler *ptr = fdlr.ptr;
  for (i = 0; i < lelt; i++) {
    ptr[i].seq = i;
    ptr[i].proc = gsc->id;
    ptr[i].part0 = 0;
  }

  int iter = 0;
  if (algo == 0)
    iter = lanczos(fiedler, &wl, initv, gsc, max_iter, &fdlr, nelg, buf);
  else if (algo == 1)
    iter = inverse(fiedler, &wl, &d, initv, gsc, max_iter, 0 /* grammian */,
                   &fdlr, nelg, buf, gid);
  metric_acc(FIEDLER_NITER, iter);

  GenmapScalar norm = 0;
  for (i = 0; i < lelt; i++)
    norm += fiedler->data[i] * fiedler->data[i];
  GenmapScalar normi;
  comm_allreduce(gsc, gs_double, gs_add, &norm, 1, &normi);
  normi = 1.0 / sqrt(norm);

  vec_scale(fiedler, fiedler, normi);

  for (i = 0; i < lelt; i++)
    elems[i].fiedler = fiedler->data[i];

  array_free(&fdlr);

  vec_destroy(initv);
  vec_destroy(fiedler);
  laplacian_free(&wl);

#if defined(GENMAP_OCCA)
  occa_lanczos_free();
#endif
  if (algo == 1) {
    mg_free(&d);
    laplacian_free(&l);
  }

  return 0;
}
