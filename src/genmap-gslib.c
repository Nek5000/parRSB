#include "genmap-impl.h"
#include "genmap-gslib.h"

int GenmapCrystalInit(GenmapHandle h, GenmapComm c) {
  crystal_init(&(h->cr), &(c->gsComm));
  return 0;
}

int GenmapCrystalTransfer(GenmapHandle h, int field) {
  if(field == GENMAP_ORIGIN)
    sarray_transfer(struct GenmapElement_private, &(h->elementArray), origin, 0,
                    &(h->cr));
  else if(field == GENMAP_PROC)
    sarray_transfer(struct GenmapElement_private, &(h->elementArray), proc, 0,
                    &(h->cr));
  return 0;
}

int GenmapCrystalFinalize(GenmapHandle h) {
  crystal_free(&(h->cr));
  return 0;
}

void GenmapSortLocal(GenmapHandle h, int field) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  if(field == GENMAP_FIEDLER) { // Fiedler
    sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt, fiedler,
                  TYPE_DOUBLE, globalId, TYPE_LONG, &h->buf);
  } else if(field == GENMAP_GLOBALID) {
    sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt,
                  globalId, TYPE_LONG, fiedler, TYPE_DOUBLE, &h->buf);
  }
}

