#include "genmap-impl.h"

int genmap_comm_size(genmap_comm c) { return (int)c->np; }

int genmap_comm_rank(genmap_comm c) { return (int)c->id; }

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
