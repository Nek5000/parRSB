#include <ctype.h>
#include <genmap-impl.h>

#include <occa.h>

#define BLK_SIZE 512

static occaDevice device;

int occa_init(char *backend_, int device_id, int platform_id) {
  size_t len = strnlen(backend_, BUFSIZ);

  char backend[BUFSIZ] = {'\0'};
  int i;
  for (i = 0; i < len; i++)
    backend[i] = tolower(backend_[i]);

  char *fmt_ocl = "{{mode: '%s'}, {device_id: %d}, {platform_id: %d}}";
  char *fmt_cuda = "{{mode: '%s'}, {device_id: %d}}";

  char fmt[BUFSIZ] = {'\0'};
  if (strncmp(backend, "opencl", BUFSIZ) == 0)
    snprintf(fmt, BUFSIZ, fmt_ocl, "OpenCL", device_id, platform_id);
  else if (strncmp(backend, "cuda", BUFSIZ) == 0)
    snprintf(fmt, BUFSIZ, fmt_cuda, "CUDA", device_id);
  else
    return 1;

  device = occaCreateDeviceFromString(fmt);

  return 0;
}

static occaKernel copy, sum, scale, dot, norm, add2s1, add2s2, addc, lplcn;
static occaMemory o_p, o_w, o_r, o_rr, o_v, o_off, o_x, o_y;
static occaJson props;

static occaMemory o_wrk;
static GenmapScalar *wrk;
static size_t wrk_size;

int occa_lanczos_init(struct comm *c, struct laplacian *l, int niter) {
  int type = l->type;
  assert((type & CSR) == CSR);

  /* device alloc */
  uint lelt = l->nel;
  o_p =
      occaDeviceMalloc(device, sizeof(GenmapScalar) * lelt, NULL, occaDefault);
  o_w =
      occaDeviceMalloc(device, sizeof(GenmapScalar) * lelt, NULL, occaDefault);
  o_r =
      occaDeviceMalloc(device, sizeof(GenmapScalar) * lelt, NULL, occaDefault);
  o_y =
      occaDeviceMalloc(device, sizeof(GenmapScalar) * lelt, NULL, occaDefault);
  o_rr = occaDeviceMalloc(device, sizeof(GenmapScalar) * lelt * (niter + 1),
                          NULL, occaDefault);

  uint wrk_size = (lelt + BLK_SIZE - 1) / BLK_SIZE;
  o_wrk = occaDeviceMalloc(device, sizeof(GenmapScalar) * wrk_size, NULL,
                           occaDefault);
  wrk = tcalloc(GenmapScalar, (niter + 1) * lelt);

  /* laplacian */
  struct csr_mat *M = l->M;
  assert(M->rn == lelt);
  uint nnz = M->roff[M->rn];

  o_v = occaDeviceMalloc(device, sizeof(GenmapScalar) * nnz, NULL, occaDefault);
  occaCopyPtrToMem(o_v, M->v, occaAllBytes, 0, occaDefault);

  o_off =
      occaDeviceMalloc(device, sizeof(uint) * (lelt + 1), NULL, occaDefault);
  occaCopyPtrToMem(o_off, M->roff, occaAllBytes, 0, occaDefault);

  o_x = occaDeviceMalloc(device, sizeof(GenmapScalar) * nnz, NULL, occaDefault);

  /* occa kernels */
  char okl[PATH_MAX];
  char *okl_dir = getenv("PARRSB_OKL_DIR");
  if (okl_dir != NULL)
    strncpy(okl, okl_dir, PATH_MAX);
  strncat(okl, "/occa.okl", PATH_MAX);

  props = occaCreateJson();
  occaJsonObjectSet(props, "defines/BLK_SIZE", occaInt(512));
  occaJsonObjectSet(props, "defines/uint", occaString("unsigned int"));
  occaJsonObjectSet(props, "defines/ulong", occaString("unsigned long"));
  occaJsonObjectSet(props, "defines/scalar", occaString("double"));

  if (c->id == 0) {
    copy = occaDeviceBuildKernel(device, okl, "copy", props);
    sum = occaDeviceBuildKernel(device, okl, "sum", props);
    scale = occaDeviceBuildKernel(device, okl, "scale", props);
    dot = occaDeviceBuildKernel(device, okl, "dot", props);
    norm = occaDeviceBuildKernel(device, okl, "norm", props);
    add2s1 = occaDeviceBuildKernel(device, okl, "add2s1", props);
    add2s2 = occaDeviceBuildKernel(device, okl, "add2s2", props);
    addc = occaDeviceBuildKernel(device, okl, "addc", props);
    lplcn = occaDeviceBuildKernel(device, okl, "laplacian_csr", props);
  }

  comm_barrier(c);

  copy = occaDeviceBuildKernel(device, okl, "copy", props);
  sum = occaDeviceBuildKernel(device, okl, "sum", props);
  scale = occaDeviceBuildKernel(device, okl, "scale", props);
  dot = occaDeviceBuildKernel(device, okl, "dot", props);
  norm = occaDeviceBuildKernel(device, okl, "norm", props);
  add2s1 = occaDeviceBuildKernel(device, okl, "add2s1", props);
  add2s2 = occaDeviceBuildKernel(device, okl, "add2s2", props);
  addc = occaDeviceBuildKernel(device, okl, "addc", props);
  lplcn = occaDeviceBuildKernel(device, okl, "laplacian_csr", props);

  return 0;
}

int occa_lanczos_aux(genmap_vector diag, genmap_vector upper, genmap_vector *rr,
                     uint lelt, ulong nelg, int niter, genmap_vector f,
                     struct laplacian *gl, struct comm *gsc, buffer *bfr) {
  occaCopyPtrToMem(o_r, f->data, occaAllBytes, 0, occaDefault);

  // vec_ortho(gsc, r, nelg);
  /* orthogonalize */
  occaKernelRun(sum, o_wrk, o_r, occaInt(lelt));
  occaCopyMemToPtr(wrk, o_wrk, occaAllBytes, 0, occaDefault);

  GenmapScalar tmp = 0.0;
  uint i;
  for (i = 0; i < wrk_size; i++)
    tmp += wrk[i];
  GenmapScalar rtr;
  comm_allreduce(gsc, gs_double, gs_add, &tmp, 1, &rtr);
  tmp /= nelg;

  occaKernelRun(addc, o_r, occaDouble(-1.0 * tmp), occaInt(lelt));

  /* normalize */
  occaKernelRun(norm, o_wrk, occaInt(lelt), o_r);
  occaCopyMemToPtr(wrk, o_wrk, occaAllBytes, 0, occaDefault);

  rtr = 0;
  for (i = 0; i < wrk_size; i++)
    rtr += wrk[i];
  comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, &tmp);
  GenmapScalar rnorm = sqrt(rtr);
  GenmapScalar rni = 1.0 / rnorm;

  // vec_scale(rr[0], r, rni);
  occaKernelRun(scale, o_rr, occaUInt(0), o_r, occaDouble(rni), occaUInt(lelt));

  GenmapScalar eps = 1.e-5;
  GenmapScalar rtol = rnorm * eps;

  GenmapScalar rtz1 = 1.0;
  GenmapScalar pap = 0.0;
  GenmapScalar alpha, beta;
  GenmapScalar rtz2, pap_old;

  uint indx;
  int iter;
  for (iter = 0; iter < niter; iter++) {
    rtz2 = rtz1;
    rtz1 = rtr;
    beta = rtz1 / rtz2;
    if (iter == 0)
      beta = 0.0;

    // add2s1(p, r, beta, n)
    occaKernelRun(add2s1, o_p, o_r, occaDouble(beta), occaInt(lelt));

    // vec_ortho(gsc, p, nelg);
    /* orthogonalize */
    occaKernelRun(sum, o_wrk, o_p, occaInt(lelt));
    occaCopyMemToPtr(wrk, o_wrk, occaAllBytes, 0, occaDefault);

    tmp = 0.0;
    for (int i = 0; i < wrk_size; i++)
      tmp += wrk[i];
    comm_allreduce(gsc, gs_double, gs_add, &tmp, 1, &rtr);
    tmp /= nelg;

    occaKernelRun(addc, o_r, occaDouble(-1.0 * tmp), occaInt(lelt));

    // laplacian
    occaCopyMemToPtr(wrk, o_p, occaAllBytes, 0, occaDefault);
    csr_mat_gather(gl->M->buf, gl->M, wrk, bfr);
    occaCopyPtrToMem(o_x, gl->M->buf, occaAllBytes, 0, occaDefault);
    occaKernelRun(lplcn, o_w, occaInt(gl->M->rn), o_off, o_v, o_x);

    pap_old = pap;

    // pap = vec_dot(w, p);
    occaKernelRun(dot, o_wrk, occaInt(lelt), o_w, o_p);
    occaCopyMemToPtr(wrk, o_wrk, occaAllBytes, 0, occaDefault);

    pap = 0.0;
    for (int i = 0; i < wrk_size; i++)
      pap += wrk[i];
    comm_allreduce(gsc, gs_double, gs_add, &pap, 1, &tmp);

#if 0
    if (gsc->id == 0)
      printf("iter = %d beta = %lf pp = %lf pap = %lf ww = %lf\n", iter,
             beta, pp, pap, ww);
#endif

    alpha = rtz1 / pap;

    // add2s2(r, w, -1.0 * alpha, n);
    occaKernelRun(add2s2, o_r, o_w, occaDouble(-1.0 * alpha), occaInt(lelt));

    // rtr = vec_dot(r, r);
    occaKernelRun(norm, o_wrk, occaInt(lelt), o_r);
    occaCopyMemToPtr(wrk, o_wrk, occaAllBytes, 0, occaDefault);

    rtr = 0;
    for (int i = 0; i < wrk_size; i++)
      rtr += wrk[i];
    comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, &tmp);
    rnorm = sqrt(rtr);
    rni = 1.0 / rnorm;

    // vec_scale(rr[iter + 1], r, rni);
    indx = (iter + 1) * lelt;
    occaKernelRun(scale, o_rr, occaUInt(indx), o_r, occaDouble(rni),
                  occaUInt(lelt));

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

  occaCopyMemToPtr(wrk, o_rr, occaAllBytes, 0, occaDefault);
  for (i = 0; i < iter + 1; i++) {
    indx = i * lelt;
    memcpy(rr[i]->data, &wrk[indx], sizeof(GenmapScalar) * lelt);
  }

  metric_acc(TOL_FINAL, rnorm);
  metric_acc(TOL_TARGET, rtol);

  return iter;
}

int occa_lanczos_free() {
  occaFree(&o_p);
  occaFree(&o_w);
  occaFree(&o_r);
  occaFree(&o_y);
  occaFree(&o_rr);

  occaFree(&o_v);
  occaFree(&o_off);
  occaFree(&o_x);

  occaFree(&o_wrk);
  if (wrk != NULL)
    free(wrk);

  occaFree(&copy);
  occaFree(&sum);
  occaFree(&scale);
  occaFree(&dot);
  occaFree(&norm);
  occaFree(&add2s1);
  occaFree(&add2s2);
  occaFree(&addc);

  return 0;
}

int occa_free() {
  occaFree(&device);
  return 0;
}
