#include <ctype.h>
#include <genmap-impl.h>

#include <occa.h>

static occaDevice device;

int occa_init(char *backend, int device_id, int platform_id) {
  size_t len = strnlen(backend, BUFSIZ);

  char be[BUFSIZ];
  int i;
  for (i = 0; i < len; i++)
    be[i] = tolower(backend[i]);

  char *fmt_ocl = "mode: %s, device_id: %d, platform_id: %d";
  char *fmt_cuda = "mode: %s, device_id: %d";

  char fmt[BUFSIZ];
  if (strncmp(backend, "opencl", BUFSIZ) == 0)
    snprintf(fmt, BUFSIZ, fmt_ocl, "OpenCL", device_id, platform_id);
  else if (strncmp(backend, "cuda", BUFSIZ) == 0)
    snprintf(fmt, BUFSIZ, fmt_cuda, "CUDA", device_id);
  else
    return 1;

  device = occaCreateDeviceFromString(fmt);
}

static occaKernel copy, sum, scale, dot, add2s1, add2s2;
static occaKernel laplacian_csr;
static occaMemory o_p, o_w, o_r, o_rr;

int occa_lanczos_setup(struct laplacian *gl, uint lelt, int niter) {
  o_r = occaDeviceMalloc(device, sizeof(scalar) * lelt, NULL, occaDefault);
  o_p = occaDeviceMalloc(device, sizeof(scalar) * lelt, NULL, occaDefault);
  o_w = occaDeviceMalloc(device, sizeof(scalar) * lelt, NULL, occaDefault);
  o_w = occaDeviceMalloc(device, sizeof(scalar) * lelt * (niter + 1), NULL, occaDefault);
}

int occa_lanczos_free() {}
