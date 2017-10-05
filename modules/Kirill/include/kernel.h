//
// Created by kirill on 23.10.16.
//

#ifndef HEAT_EQUATION_SPARSE_SPARSE_H
#define HEAT_EQUATION_SPARSE_SPARSE_H

#include "ts.h"

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

// Суммирование векторов для рунге-кутты
void sumV(TYPE **result, TYPE *U, TYPE *k1, TYPE *k2, TYPE *k3, TYPE *k4, int N, TYPE h);

void printSpMat(SpMatrix mat);
TYPE procedure(SpMatrix mat, int i, int j);

// Умножение плотной матрицы на вектор
void denseMult(TYPE **result, TYPE **mat, TYPE *vec, int dim);

// Создание матрицы для явных схем
void createExplicitSpMat(SpMatrix *mat, TYPE coeffs[4], int dim, int NX, int NXY);

void createExplicitSpMatV2R(SpMatrix* mat, TYPE* coeffs, int nx, int ny, int nz, MPI_Comm comm);
void createExplicitSpMatV2(SpMatrix *mat, TYPE coeffs[4], int nx, int ny, int nz);

// Создание матрицы для неявных схем
void createImplicitSpMat(SpMatrix* mat, TYPE coeffs[3], int nx, int ny, int nz);

#ifdef __cplusplus
}
#endif


#endif //HEAT_EQUATION_SPARSE_SPARSE_H
