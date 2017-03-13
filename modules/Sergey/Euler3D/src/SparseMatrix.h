//
// Created by lenferd on 27.10.16.
//

#ifndef SPARSEMATRIX_SPARSEMATRIX_H
#define SPARSEMATRIX_SPARSEMATRIX_H
#include <omp.h>
#include <cstdio>

const int ENABLE_PARALLEL = 1;

typedef struct SparseMatrix {
    int _size;
    int _rows;
    double *values;
    int *columns;
    int *pointerB;
} SpaceMatrix;

void fillMatrix2Expr(SparseMatrix &sp, int size, double expr1, double expr2);
void multiplicateVector(SparseMatrix &sp, double *&vect, double *&result, int size);
void spMatrixInit(SparseMatrix &sp, int size, int rows);
void printVectors(SparseMatrix &sp);

#endif //SPARSEMATRIX_SPARSEMATRIX_H
