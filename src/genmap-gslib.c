#include "genmap-impl.h"
#include "genmap-gslib.h"
//
// comm, buffer and array init, finalize
//
int GenmapGSCommInit(GenmapComm c, GenmapCommExternal ce) {
  comm_init(&c->gsComm, ce);
  c->verticesHandle =  NULL;
}

int GenmapGSCommFinalize(GenmapComm c) {
  if(&c->gsComm)
    comm_free(&c->gsComm);
}

int GenmapBufferInit(GenmapComm c, size_t size) {
  buffer_init(&c->buf, size);
}

int GenmapBufferFinalize(GenmapComm c) {
  if(&c->buf)
    buffer_free(&c->buf);
}
// Crystal Router
//
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
//
// Sort functionality
//
void GenmapSortLocal(GenmapHandle h, int field) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  // To get the buffer space
  GenmapComm c = GenmapGetGlobalComm(h);

  if(field == GENMAP_FIEDLER) { // Fiedler
    sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt,
                  fiedler, TYPE_DOUBLE, globalId, TYPE_LONG, &c->buf);
  } else if(field == GENMAP_GLOBALID) {
    sarray_sort_2(struct GenmapElement_private, elements, (GenmapUInt)lelt,
                  globalId, TYPE_LONG, fiedler, TYPE_DOUBLE, &c->buf);
  }
}
//
// gather-scatter functionality
//
int GenmapGSSetup(GenmapComm c, GenmapLong *vertices, GenmapUInt numPoints) {
  if(c->verticesHandle)
    gs_free(c->verticesHandle);

  c->verticesHandle = gs_setup(vertices, numPoints, &c->gsComm, 0,
                               gs_crystal_router, 0);
  return c->verticesHandle != NULL;
}

void GenmapGSAdd(GenmapComm c, GenmapScalar *u) {
  gs(u, genmap_gs_scalar, gs_add, 0, c->verticesHandle, &c->buf);
}

void GenmapGSMax(GenmapComm c, GenmapLong *u) {
  gs(u, genmap_gs_long, gs_max, 0, c->verticesHandle, &c->buf);
}

int GenmapGSFree(GenmapComm c) {
  if(c->verticesHandle)
    gs_free(c->verticesHandle);
  return 0;
}
//
// Scan functionality
//
void GenmapScan(GenmapHandle h, GenmapComm c) {
  GenmapLong out[2][1], buf[2][1];
  GenmapLong lelt = GenmapGetNLocalElements(h);
  comm_scan(out, &(c->gsComm), genmap_gs_long, gs_add, &lelt, 1, buf);
  GenmapSetLocalStartIndex(h, out[0][0]);
  GenmapSetNGlobalElements(h, out[1][0]);
}
