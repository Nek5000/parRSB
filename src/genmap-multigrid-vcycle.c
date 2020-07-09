#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

void mg_vcycle(GenmapScalar *u,GenmapScalar *rhs,mgData d){

  uint lvl;
  for(lvl=0; lvl<d->nLevels-1; lvl++){
    parMat M=d->levels[lvl]->M;
    uint off=d->level_off[lvl];
    uint n  =d->level_off[lvl+1]-off;
  }
}
