//
// Created by lenferd on 27.10.16.
//

#ifndef SPARSEMATRIX_SPARSEMATRIX_H
#define SPARSEMATRIX_SPARSEMATRIX_H
#include <omp.h>
#include <cstdio>
#include "Task.h"
#include "StructDeclamer.h"

const int ENABLE_PARALLEL = 1;

struct MatrixValue;

struct SparseMatrix {
    int _size;
    int _rows;
    double *values;
    int *columns;   // какой столбец
    int *pointerB;  // указатель на начало строки
};


void fillMatrix2Expr(SparseMatrix &sp, int size, double expr1, double expr2);

void fillMatrix3d6Expr(SparseMatrix &sp, MatrixValue &taskexpr, int sizeX, int sizeY, int sizeZ);
void fillMatrix3d6Expr_wo_boundaries(SparseMatrix &sp, MatrixValue &taskexpr, int sizeX, int sizeY, int sizeZ);

void multiplicateVector(SparseMatrix &sp, double *&vect, double *&result, int size);
void spMatrixInit(SparseMatrix &sp, int size, int rows);
void printVectors(SparseMatrix &sp);


#endif //SPARSEMATRIX_SPARSEMATRIX_H
