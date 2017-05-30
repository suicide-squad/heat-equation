//
// Created by kirill on 23.10.16.
//

#ifndef HEAT_EQUATION_SPARSE_SPARSE_H
#define HEAT_EQUATION_SPARSE_SPARSE_H

#define ENABLE_PARALLEL 0
//#define _COMPLEX_

#include <omp.h>
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif


#ifdef _COMPLEX_

#include <complex.h>
typedef complex TYPE;

#else
typedef double TYPE;
#endif

// CSR (Compressed Sparse Rows)
typedef struct {
  double* value;    // Элементы матрицы
  int* col;         // Номера столбцов для каждого элемента
  int* rowIndex;    // Место каждого ненулевого элемента в каждой строке
  size_t nz;        // Количество ненулевых элементов
  size_t nRows;     // Количество строк
} SpMatrix;

void initSpMat(SpMatrix* mat, size_t nz, size_t nRows);
void freeSpMat(SpMatrix* mat);

// Умножение разреженной матрицы на вектор
void multMV(TYPE** result, SpMatrix matrix, TYPE* vector);

// Суммирование векторов для рунге-кутты
void sumV(TYPE **result, TYPE *U, TYPE *k1, TYPE *k2, TYPE *k3, TYPE *k4, size_t N, double h);

void printSpMat(SpMatrix mat);
TYPE procedure(SpMatrix mat, int i, int j);

// Умножение плотной матрицы на вектор
void denseMult(TYPE **result, TYPE **mat, TYPE *vec, size_t dim);

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
