#include "genmap.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _PARRSB_OCCA_H_
#define _PARRSB_OCCA_H_

int parRSBOccaSetup(GenmapHandle h);
int parRSBOccaLaplacianSetup(GenmapHandle h);
int parRSBOccaLaplacian(GenmapHandle h);
int parRSBOccaLanczosSetup(GenmapHandle h);
int parRSBOccaLanczos(GenmapHandle h);

#endif

#ifdef __cplusplus
}
#endif
