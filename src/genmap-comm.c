#include "genmap-impl.h"

int GenmapCreateComm(genmap_comm *c_, comm_ext ce) {
  GenmapMalloc(1, c_);
  genmap_comm c = *c_;
  comm_init(&c->gsc, ce);

  c->gsh = NULL;
  c->M = NULL;

  return 0;
}

int GenmapDestroyComm(genmap_comm c) {
  if (c->gsh)
    gs_free(c->gsh);
  if (c->M)
    csr_mat_free(c->M);

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

int GenmapGop(genmap_comm c, void *v, GenmapInt size, GenmapDataType type,
              GenmapInt op) {
  if (op == GENMAP_SUM) {
    MPI_Allreduce(MPI_IN_PLACE, v, size, type, MPI_SUM, c->gsc.c);
  } else if (op == GENMAP_MAX) {
    MPI_Allreduce(MPI_IN_PLACE, v, size, type, MPI_MAX, c->gsc.c);
  } else if (op == GENMAP_MIN) {
    MPI_Allreduce(MPI_IN_PLACE, v, size, type, MPI_MIN, c->gsc.c);
  }
  return 0;
}

int GenmapReduce(genmap_comm c, void *out, void *in, GenmapInt size,
                 GenmapDataType type, GenmapInt op) {
  if (op == GENMAP_SUM) {
    MPI_Reduce(in, out, size, type, MPI_SUM, 0, c->gsc.c);
  } else if (op == GENMAP_MAX) {
    MPI_Reduce(in, out, size, type, MPI_MAX, 0, c->gsc.c);
  } else if (op == GENMAP_MIN) {
    MPI_Reduce(in, out, size, type, MPI_MIN, 0, c->gsc.c);
  }
  return 0;
}

int GenmapBcast(genmap_comm c, void *in, GenmapInt count, GenmapDataType type) {
  return MPI_Bcast(in, count, type, 0, c->gsc.c);
}

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

void GenmapSplitComm(genmap_handle h, genmap_comm *c, int bin) {
  comm_ext local;
  int id = genmap_comm_rank(*c);
  MPI_Comm_split((*c)->gsc.c, bin, id, &local);
  GenmapCrystalFinalize(h);
  GenmapDestroyComm(*c);
  GenmapCreateComm(c, local);
  MPI_Comm_free(&local);
  GenmapCrystalInit(h, *c);
}

int GenmapCrystalInit(genmap_handle h, genmap_comm c) {
  crystal_init(&(h->cr), &(c->gsc));
  return 0;
}

int GenmapCrystalTransfer(genmap_handle h, int field) {
  if (field == GENMAP_ORIGIN)
    sarray_transfer(struct rsb_element, h->elements, origin, 0, &h->cr);
  else if (field == GENMAP_PROC)
    sarray_transfer(struct rsb_element, h->elements, proc, 0, &h->cr);
  return 0;
}

int GenmapCrystalFinalize(genmap_handle h) {
  crystal_free(&(h->cr));
  return 0;
}

void comm_split(struct comm *old, int bin, int key, struct comm *new) {
#ifdef MPI
  MPI_Comm new_comm;
  MPI_Comm_split(old->c, bin, key, &new_comm);
  comm_init(new, new_comm);
  MPI_Comm_free(&new_comm);
#else
  comm_init(new, 1);
#endif
}
