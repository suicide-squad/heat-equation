//
// Created by lenferd on 27.10.16.
//

#include "SparseMatrix.h"


void spMatrixInit(SparseMatrix &sp, int size, int rows) {
    sp._size = size;
    sp._rows = rows;
    sp.values = new double[size];
    sp.columns = new int[size];
    sp.pointerB = new int[rows+1];
}

void multiplicateVector(SparseMatrix &sp, double *&vect, double *&result, int size) {

    omp_set_num_threads(4);

    #pragma omp parallel for if (ENABLE_PARALLEL)
    for (int i = 0; i < size; i++){  // iteration FOR RESULT VECTOR!!!
        double local_result = 0;
        for (int j = sp.pointerB[i]; j < sp.pointerB[i+1]; j++) {
            local_result += sp.values[j] * vect[sp.columns[j]];
        }
        result[i] = local_result;
    }
}

void fillMatrix2Expr(SparseMatrix &sp, int size, double expr1, double expr2) {
    int index = 0;
    int pIndex = 0;

    sp.values[index] = 1;
    sp.columns[index] = 0;
    sp.pointerB[pIndex++] = 0;
    ++index;

    for (int i = 1; i < size - 1; ++i) {

        //printf("index %d \n", index);
        sp.values[index] = expr1;
        sp.columns[index] = i - 1;
        sp.pointerB[pIndex++] = index;
        ++index;

        sp.values[index] = expr2;
        sp.columns[index] = i;
        ++index;

        sp.values[index] = expr1;
        sp.columns[index] = i + 1;
        ++index;
    }

    sp.values[index] = 1;
    sp.columns[index] = size - 1;
    sp.pointerB[pIndex++] = index;

    sp.pointerB[pIndex] = index + 1;   //end
}


void printVectors(SparseMatrix &sp) {
    printf("values\n");
    for (int i = 0; i < sp._size; ++i) {
        printf("%lf ", sp.values[i]);
    }
    printf("\n");

    printf("columns\n");
    for (int i = 0; i < sp._size; ++i) {
        printf("%d ", sp.columns[i]);
    }
    printf("\n");

    printf("pointerB\n");
    for (int i = 0; i < sp._rows + 1; ++i) {
        printf("%d ", sp.pointerB[i]);
    }
    printf("\n");
}
