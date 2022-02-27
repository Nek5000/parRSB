#include "coarse.h"
#include "genmap-impl.h"
#include "multigrid.h"
#include "sort.h"
#include <math.h>
#include <time.h>

#define MM 500

static int vec_create(genmap_vector *x, GenmapInt size) {
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

inline static scalar dot(scalar *y, scalar *x, uint n) {
  scalar result = 0.0;
  for (uint i = 0; i < n; i++)
    result += x[i] * y[i];

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
                   struct mg_data *d, struct comm *c, int miter, int null_space,
                   int verbose, buffer *bfr) {
  slong out[2][1], buf[2][1], in = n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
  ulong ng = out[1][0];

  if (ng == 0)
    return 0;

  scalar *z = (scalar *)tcalloc(scalar, (6 + 2 * miter + 2) * n);
  scalar *w = z + n, *r = w + n, *p = r + n, *z0 = p + n, *dz = z0 + n;
  scalar *P = dz + n, *W = P + n * (miter + 1);

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
  struct mg_data *d = mg_setup(L, factor, gsc, &cr, buf);
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

    // add2s1(p,r,beta,n)
    uint i;
    for (i = 0; i < lelt; i++)
      p->data[i] = beta * p->data[i] + r->data[i];

    GenmapScalar pp = vec_dot(p, p);
    comm_allreduce(gsc, gs_double, gs_add, &pp, 1, buf);

    vec_ortho(gsc, p, nelg);

    metric_tic(gsc, LAPLACIAN);
    laplacian(w->data, gl, p->data, bfr);
    metric_toc(gsc, LAPLACIAN);

    GenmapScalar ww = vec_dot(w, w);
    comm_allreduce(gsc, gs_double, gs_add, &ww, 1, buf);

    pap_old = pap, pap = vec_dot(w, p);
    comm_allreduce(gsc, gs_double, gs_add, &pap, 1, buf);
#if 0
    if (gsc->id == 0)
      printf("host iter = %d beta = %lf pp = %lf pap = %lf\n", iter, beta, pp,
             pap);
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
      upper->size = iter++;
      diag->size = iter;
      break;
    }
  }

  metric_acc(TOL_FNL, rnorm);
  metric_acc(TOL_TGT, rtol);

  vec_destroy(p);
  vec_destroy(w);
  vec_destroy(r);

  return iter;
}

static int lanczos(genmap_vector fiedler, struct rsb_element *elems, int nv,
                   genmap_vector initv, struct comm *gsc, int max_iter,
                   slong nelg, buffer *bfr) {
  assert(fiedler->size == initv->size);
  uint lelt = fiedler->size;

  struct laplacian *wl;
#if defined(GENMAP_OCCA)
  wl = laplacian_init(elems, lelt, nv, CSR, gsc, bfr);
  occa_lanczos_init(gsc, wl, max_iter);
#else
  wl = laplacian_init(elems, lelt, nv, CSR, gsc, bfr);
#endif

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

  uint i;
  for (i = 0; i < max_iter + 1; ++i)
    vec_create(&rr[i], lelt);

  int iter, ipass = 0;
  do {
#if defined(GENMAP_OCCA)
    iter = occa_lanczos_aux(alpha, beta, rr, lelt, nelg, max_iter, initv, wl,
                            gsc, bfr);
#else
    iter =
        lanczos_aux(alpha, beta, rr, lelt, nelg, max_iter, initv, wl, gsc, bfr);
#endif

    genmap_vector evTriDiag;
    vec_create(&evTriDiag, iter);

    // Use TQLI and find the minimum eigenvalue and associated vector
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
#if 0
    if (gsc->id == 0)
      printf("ipass = %d eValMinI = %d eValue = %lf\n", ipass, eValMinI,
             eValues->data[eValMinI]);
#endif

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

#if defined(GENMAP_OCCA)
  occa_lanczos_free();
#endif

  laplacian_free(wl);

  return ipass;
}

int fiedler(struct rsb_element *elems, uint lelt, int nv, int max_iter,
            int algo, struct comm *gsc, buffer *buf) {
  slong out[2][1], bfr[2][1], in = lelt;
  comm_scan(out, gsc, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0], nelg = out[1][0];

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

  int iter = 0;
  if (algo == 0)
    iter = lanczos(fiedler, elems, nv, initv, gsc, max_iter, nelg, buf);
  else if (algo == 1)
    iter = inverse(fiedler->data, fiedler->size, elems, nv, initv->data, gsc,
                   max_iter, nelg, buf);
  metric_acc(FIEDLER_NITER, iter);

  GenmapScalar norm = 0;
  for (i = 0; i < lelt; i++)
    norm += fiedler->data[i] * fiedler->data[i];
#if 0
  if (gsc->id == 0)
    printf("fiedler norm = %lf\n", norm);
#endif

  GenmapScalar normi;
  comm_allreduce(gsc, gs_double, gs_add, &norm, 1, &normi);
  normi = 1.0 / sqrt(norm);

  vec_scale(fiedler, fiedler, normi);

  for (i = 0; i < lelt; i++)
    elems[i].fiedler = fiedler->data[i];

  vec_destroy(initv);
  vec_destroy(fiedler);

  return 0;
}
