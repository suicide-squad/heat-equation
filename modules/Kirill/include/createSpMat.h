//
// Created by kirill on 12.11.17.
//

#ifndef HEAT_EQUATION_CREATESPMAT_H
#define HEAT_EQUATION_CREATESPMAT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "utils/ts.h"
#include "multMV.h"

// Создание матрицы для явных схем
void createExplicitSpMat(SpMatrix *mat, TYPE coeffs[4], int dim, int NX, int NXY);

void createExplicitSpMatV2R(SpMatrix* mat, const TYPE* const coeffs, const int nx, const int ny, const int nz,
                            const int rank_left, const int rank_right, const int rank_down, const int rank_top);
void createExplicitSpMatV2(SpMatrix *mat, TYPE coeffs[4], int nx, int ny, int nz);

// Создание матрицы для неявных схем
void createImplicitSpMat(SpMatrix* mat, TYPE coeffs[3], int nx, int ny, int nz);

#ifdef __cplusplus
}
#endif

#endif //HEAT_EQUATION_CREATESPMAT_H
