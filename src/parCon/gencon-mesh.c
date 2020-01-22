#include <stdlib.h>

#include "gencon-impl.h"

int freeMesh(Mesh mesh){
  exaArrayFree(mesh->elements);
  free(mesh);
  mesh=NULL;
}
