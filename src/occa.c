#include <ctype.h>
#include <genmap-impl.h>

#include <occa.h>

#define scalar GenmapScalar

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

static occaKernel copy, sum, scale, dot, add2s1, add2s2;
static occaKernel laplacian_csr;
static occaMemory o_p, o_w, o_r, o_rr;
static occaJson props;

int occa_lanczos_init(struct laplacian *l, int niter) {
  int type = l->type;
  assert((type & CSR) == CSR);

  uint lelt = l->nel;
  o_r = occaDeviceMalloc(device, sizeof(scalar) * lelt, NULL, occaDefault);
  o_p = occaDeviceMalloc(device, sizeof(scalar) * lelt, NULL, occaDefault);
  o_w = occaDeviceMalloc(device, sizeof(scalar) * lelt, NULL, occaDefault);
  o_rr = occaDeviceMalloc(device, sizeof(scalar) * lelt * (niter + 1), NULL,
                          occaDefault);

  props = occaCreateJson();
  copy = occaDeviceBuildKernel(device, "occa.okl", "copy", props);
  sum = occaDeviceBuildKernel(device, "occa.okl", "sum", props);
  scale = occaDeviceBuildKernel(device, "occa.okl", "scale", props);
  dot = occaDeviceBuildKernel(device, "occa.okl", "dot", props);
  add2s1 = occaDeviceBuildKernel(device, "occa.okl", "add2s1", props);
  add2s2 = occaDeviceBuildKernel(device, "occa.okl", "add2s2", props);

  return 0;
}

int occa_lanczos(scalar *diag, scalar *upper, scalar *rr, scalar *f,
                 struct comm *c, buffer *bfr) {
  return 0;
}

int occa_lanczos_free() {
  occaFree(&o_p);
  occaFree(&o_w);
  occaFree(&o_r);
  occaFree(&o_rr);
  occaFree(&device);
  return 0;
}
