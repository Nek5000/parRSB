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
static occaMemory o_p, o_w, o_r, o_rr;
static occaJson props;

int occa_lanczos_init(struct comm *c, struct laplacian *l, int niter) {
  int type = l->type;
  assert((type & CSR) == CSR);

  uint lelt = l->nel;
  o_r =
      occaDeviceMalloc(device, sizeof(GenmapScalar) * lelt, NULL, occaDefault);
  o_p =
      occaDeviceMalloc(device, sizeof(GenmapScalar) * lelt, NULL, occaDefault);
  o_w =
      occaDeviceMalloc(device, sizeof(GenmapScalar) * lelt, NULL, occaDefault);
  o_rr = occaDeviceMalloc(device, sizeof(GenmapScalar) * lelt * (niter + 1),
                          NULL, occaDefault);

  char okl[PATH_MAX];
  char *okl_dir = getenv("PARRSB_OKL_DIR");
  if (okl_dir != NULL)
    strncpy(okl, okl_dir, PATH_MAX);
  strncat(okl, "/occa.okl", PATH_MAX);

  props = occaCreateJson();
  if (c->id == 0) {
    copy = occaDeviceBuildKernel(device, okl, "copy", props);
    sum = occaDeviceBuildKernel(device, okl, "sum", props);
    scale = occaDeviceBuildKernel(device, okl, "scale", props);
    dot = occaDeviceBuildKernel(device, okl, "dot", props);
    add2s1 = occaDeviceBuildKernel(device, okl, "add2s1", props);
    add2s2 = occaDeviceBuildKernel(device, okl, "add2s2", props);
  }

  comm_barrier(c);

  return 0;
}

int occa_lanczos(GenmapScalar *diag, GenmapScalar *upper, GenmapScalar *rr,
                 GenmapScalar *f, struct comm *c, buffer *bfr) {
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
