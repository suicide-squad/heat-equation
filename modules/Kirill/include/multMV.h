//
// Created by kirill on 23.10.16.
//

#ifndef HEAT_EQUATION_SPARSE_SPARSE_H
#define HEAT_EQUATION_SPARSE_SPARSE_H

#include "utils/ts.h"

#include <omp.h>
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif


// CSR (Compressed Sparse Rows)
typedef struct {
  TYPE* value;    // Элементы матрицы
  int* col;         // Номера столбцов для каждого элемента
  int* rowIndex;    // Место каждого ненулевого элемента в каждой строке
  int nz;        // Количество ненулевых элементов
  int nRows;     // Количество строк
} SpMatrix;

void initSpMat(SpMatrix* mat, int nz, int nRows);
void freeSpMat(SpMatrix* mat);

void copyingBorders(TYPE* vec, int nx, int ny, int nz);

// Умножение разреженной матрицы на вектор
void multMV(TYPE* result, SpMatrix mat, TYPE* vec, int nx, int ny, int nz, TYPE* coeff);

void multMV_altera(TYPE* result, SpMatrix mat, TYPE* vec, int sizeTime);
void naive_formula(TYPE* result, TYPE* vec, const TYPE* const coeff, const int nx, const int ny, const int nz, const int sizeTime);

// Суммирование векторов для рунге-кутты
void sumV(TYPE **result, TYPE *U, TYPE *k1, TYPE *k2, TYPE *k3, TYPE *k4, int N, TYPE h);

void printSpMat(SpMatrix mat);
TYPE procedure(SpMatrix mat, int i, int j);

// Умножение плотной матрицы на вектор
void denseMult(TYPE **result, TYPE **mat, TYPE *vec, int dim);

#ifdef __cplusplus
}
#endif


#endif //HEAT_EQUATION_SPARSE_SPARSE_H
