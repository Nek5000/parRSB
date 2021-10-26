#include <ctype.h>
#include <genmap-impl.h>

#if defined(GENMAP_OCCA)

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

int occa_ax() {}

#endif
