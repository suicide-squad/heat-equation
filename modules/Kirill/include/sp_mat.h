//
// Created by kirill on 23.10.16.
//

#ifndef SPARSE_SPARSE_H
#define SPARSE_SPARSE_H

#define ENABLE_PARALLEL 1
//#define _COMPLEX_

#include <omp.h>

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
  TYPE* value;      // Элементы матрицы
  int* col;         // Номера столбцов для каждого элемента
  int* rowIndex;    // Место каждого ненулевого элемента в каждой строке
  size_t nz;        // Количество ненулевых
  size_t nRows;     // Количество строк
} spMatrix;

void initSpMat(spMatrix* mat, size_t nz, size_t nRows);
void freeSpMat(spMatrix* mat);

void multMV(TYPE** result, spMatrix matrix, TYPE* vector);


void sumV(size_t N, double h, TYPE **result, TYPE *U, TYPE *k1, TYPE *k2, TYPE *k3, TYPE *k4);

void printSpMat(spMatrix mat);
TYPE procedure(spMatrix mat, int i, int j);

#ifdef __cplusplus
}
#endif


#endif //SPARSE_SPARSE_H
