#include "genmap-impl.h"

int GenmapCreateComm(genmap_comm *c_, comm_ext ce) {
  GenmapMalloc(1, c_);
  genmap_comm c = *c_;
  comm_init(&c->gsc, ce);

  return 0;
}

int GenmapDestroyComm(genmap_comm c) {
  comm_free(&c->gsc);
  GenmapFree(c);

  return 0;
}

int genmap_comm_size(genmap_comm c) { return (int)c->gsc.np; }

int genmap_comm_rank(genmap_comm c) { return (int)c->gsc.id; }

genmap_comm GenmapGetLocalComm(genmap_handle h) { return h->local; }

void GenmapSetLocalComm(genmap_handle h, genmap_comm c) { h->local = c; }

genmap_comm GenmapGetGlobalComm(genmap_handle h) { return h->global; }

void GenmapSetGlobalComm(genmap_handle h, genmap_comm c) { h->global = c; }

void comm_split(struct comm *old, int bin, int key, struct comm *new_) {
#ifdef MPI
  MPI_Comm new_comm;
  MPI_Comm_split(old->c, bin, key, &new_comm);
  comm_init(new_, new_comm);
  MPI_Comm_free(&new_comm);
#else
  comm_init(new_, 1);
#endif
}
