#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sp_mat.h>
#include <mpi.h>

#include "sp_mat.h"

void initSpMat(SpMatrix *mat, size_t nz, size_t nRows) {
  mat->nz = nz;
  mat->nRows = nRows;
  mat->value = (TYPE *)malloc(sizeof(TYPE) * nz);
  mat->col = (int *)malloc(sizeof(int) * nz);
  mat->rowIndex = (int *)malloc(sizeof(int) * (nRows) + 1);
  memset(mat->rowIndex, 0, nRows + 1);
}

void freeSpMat(SpMatrix* mat) {
  free(mat->value);
  free(mat->col);
  free(mat->rowIndex);
}

void multMV(TYPE** result, SpMatrix mat, TYPE* vec) {
  TYPE localSum;
  #pragma omp parallel private(localSum) if (ENABLE_PARALLEL)
  {
    #pragma omp for nowait
    for (int i = 0; i < mat.nRows; i++) {
      localSum = 0.0;
      for (int j = mat.rowIndex[i]; j < mat.rowIndex[i + 1]; j++)
        localSum += mat.value[j] * vec[mat.col[j]];
      (*result)[i] = localSum;
    }
  }
}

void sumV(TYPE **result, TYPE *U, TYPE *k1, TYPE *k2, TYPE *k3, TYPE *k4, size_t N, double h) {
  #pragma omp parallel for if (ENABLE_PARALLEL)
  for (int i = 0; i < N; i++)
    (*result)[i] = U[i] + h*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
}


void printSpMat(SpMatrix mat) {
  for (int i = 0; i < mat.nRows; i++) {
    for (int j = 0; j < mat.nRows; j++)
      printf("%.0lf", procedure(mat, i, j));
    printf("\n");
  }
}

TYPE procedure(SpMatrix mat, int i, int j) {
  TYPE result = 0;
  int N1 = mat.rowIndex[i];
  int N2 = mat.rowIndex[i+1];
  for(int k = N1; k < N2; k++) {
    if (mat.col[k] == j) {
      result = mat.value[k];
      break;
    }
  }
  return result;
}

void denseMult(double **result, double **mat, double *vec, size_t dim) {
  memset(*result, 0, dim*sizeof(double));
  for (int x = 0; x < dim; x++) {
    for (int i = 0;i < dim;i++)
      (*result)[x]+=mat[x][i]*vec[i];
  }
}

void createExplicitSpMat(SpMatrix *mat, TYPE coeffs[4], int dim, int NX, int NXY) {
  int index = 0, j, k, shiftIndex;
  mat->rowIndex[0] = 0;
  for (int i = 0; i < dim; i++) {
    if (i % NX == 0)
      shiftIndex = 1;
    else if (i % NX == NX - 1)
      shiftIndex = -1;
    else
      shiftIndex = 0;

    mat->col[index] = i + shiftIndex;
    mat->value[index] = coeffs[1];
    index++;

    // ***************************************
    //      Смещение на x - 1
    // ***************************************
    mat->col[index] = i + shiftIndex - 1;
    mat->value[index] = coeffs[0];
    index++;

    // ***************************************
    //      Смещение на x + 1
    // ***************************************
    mat->col[index] = i + shiftIndex + 1;
    mat->value[index] = coeffs[0];
    index++;
    // ***************************************
    //      Смещение на y - 1
    // ***************************************
    j = (i + shiftIndex) % NXY - NX;
    if (j < 0) {
//      mat->col[index] = DIM + j;
      mat->col[index] = NXY + j;
      mat->value[index] = coeffs[2];

      index++;
    } else {
      mat->col[index] = j;
      mat->value[index] = coeffs[2];

      index++;
    }
    // ***************************************
    //      Смещение на y + 1
    // ***************************************
    mat->col[index] = ((i + shiftIndex) / NXY) * NXY + (i + shiftIndex + NX) % NXY;
    //mat->col[index] = (i + shiftIndex + NX) % NXY;
    //mat->col[index] = (i + NX) % NXY;
    mat->value[index] = coeffs[2];

    index++;
    // ***************************************
    //    Смещение на  z - 1
    // ***************************************
    k = i + shiftIndex - NXY;
    if (k <= 0) {
      mat->col[index] = dim + k;
      mat->value[index] = coeffs[3];
      index++;
    } else {
      mat->col[index] = k;
      mat->value[index] = coeffs[3];
      index++;
    }
    // ***************************************
    //    Смещение на  z + 1
    // ***************************************
    mat->col[index] = (i + shiftIndex + NXY) % dim;
    mat->value[index] = coeffs[3];
    index++;
    // ***************************************

    mat->rowIndex[i + 1] = mat->rowIndex[i] + 7;
  }
}

void createExplicitSpMatV2(SpMatrix *mat, TYPE coeffs[4], int nx, int ny, int nz) {
  int index = 0, k = 0;
  int shiftIndexX, shiftIndexY, shiftIndexZ;
  int shift;
  mat->rowIndex[0] = 0;

  int realIndex;
  for (int z = 0; z < nz; z++) {
    if (z == 0)
      shiftIndexZ = nx*ny;
    else if (z == nz - 1)
      shiftIndexZ = -nx*ny;
    else
      shiftIndexZ = 0;
    for (int y = 0; y < ny; y++) {
      if (y == 0)
        shiftIndexY = nx;
      else if (y == ny - 1)
        shiftIndexY = -nx;
      else
        shiftIndexY = 0;

      for (int x = 0; x < nx; x++) {
        if (x == 0)
          shiftIndexX = 1;
        else if (x == nx - 1)
          shiftIndexX = -1;
        else
          shiftIndexX = 0;

        shift = shiftIndexZ + shiftIndexY + shiftIndexX;
        realIndex = x + y*nx + z*nx*ny + shift;

        // ***************************************
        //      Смещение на z - 1
        // ***************************************
        mat->col[index] = realIndex - nx*ny;
        mat->value[index] = coeffs[3];
        index++;

        // ***************************************
        //      Смещение на y - 1
        // ***************************************
        mat->col[index] = realIndex - nx;
        mat->value[index] = coeffs[2];
        index++;

        // ***************************************
        //      Смещение на x - 1
        // ***************************************
        mat->col[index] = realIndex - 1;
        mat->value[index] = coeffs[1];
        index++;

// Отсутствие смещения
        mat->col[index] = realIndex;
        mat->value[index] = coeffs[0];
        index++;


        // ***************************************
        //      Смещение на x + 1
        // ***************************************
        mat->col[index] = realIndex + 1;
        mat->value[index] = coeffs[1];
        index++;


        // ***************************************
        //      Смещение на y + 1
        // ***************************************
        mat->col[index] = realIndex + nx;
        mat->value[index] = coeffs[2];
        index++;

        // ***************************************
        //      Смещение на z + 1
        // ***************************************
        mat->col[index] = realIndex + nx*ny;
        mat->value[index] = coeffs[3];
        index++;

        k++;
        mat->rowIndex[k] = mat->rowIndex[k-1] + 7;

      }

    }
  }
}

void createExplicitSpMatV2R(SpMatrix* mat, TYPE* coeffs, int nx, int ny, int nz, MPI_Comm comm) {
  int index = 0, k = 0;
  int shiftIndexX=0, shiftIndexY=0, shiftIndexZ=0;
  int shift;
  mat->rowIndex[0] = 0;

  // Определения соседних ранков в декардовой решётке
  int rank_left, rank_right, rank_down, rank_top;
  MPI_Cart_shift(comm, 1, 1, &rank_left, &rank_right);
  MPI_Cart_shift(comm, 0, 1, &rank_down, &rank_top);

  int realIndex;
  for (int z = 0; z < nz; z++) {
    if (z == 0 && rank_down==-2)
      shiftIndexZ = 4*nx*ny;
    else if (z==1&& rank_down==-2)
      shiftIndexZ = 3*nx*ny;
    else if (z==2&& rank_down==-2)
      shiftIndexZ = 2*nx*ny;
    else if (z==3&& rank_down==-2)
      shiftIndexZ = 1*nx*ny;
    else if (z == nz - 1&& rank_top==-2)
      shiftIndexZ = -1*nx*ny;
    else if (z == nz - 2&& rank_top==-2)
      shiftIndexZ = -2*nx*ny;
    else if (z == nz - 3&& rank_top==-2)
      shiftIndexZ = -3*nx*ny;
    else if (z == nz - 4&& rank_top==-2)
      shiftIndexZ = -4*nx*ny;
    else
      shiftIndexZ = 0;
    for (int y = 0; y < ny; y++) {
      if (y == 0&&rank_left==-2)
        shiftIndexY = 4*nx;
      else if (y==1&&rank_left==-2)
        shiftIndexY = 3*nx;
      else if (y==2&&rank_left==-2)
        shiftIndexY = 2*nx;
      else if (y==3&&rank_left==-2)
        shiftIndexY = 1*nx;
      else if (y == ny - 1 &&rank_right==-2)
        shiftIndexY = -1*nx;
      else if (y == ny - 2&&rank_right==-2)
        shiftIndexY = -2*nx;
      else if (y == ny - 3&&rank_right==-2)
        shiftIndexY = -3*nx;
      else if (y == ny - 4&&rank_right==-2)
        shiftIndexY = -4*nx;
      else
        shiftIndexY = 0;

      for (int x = 0; x < nx; x++) {
        if (x == 0)
          shiftIndexX = 4;
        else if (x==1)
          shiftIndexX = 3;
        else if (x==2)
          shiftIndexX = 2;
        else if (x==3)
          shiftIndexX = 1;
        else if (x == nx - 1)
          shiftIndexX = -1;
        else if (x == nx - 2)
          shiftIndexX = -2;
        else if (x == nx - 3)
          shiftIndexX = -3;
        else if (x == nx - 4)
          shiftIndexX = -4;
        else
          shiftIndexX = 0;

        shift = shiftIndexZ + shiftIndexY + shiftIndexX;
        realIndex = x + y*nx + z*nx*ny + shift;

        // ***************************************
        //      Смещение на z - 1
        // ***************************************
        mat->col[index] = realIndex - nx*ny;
        mat->value[index] = coeffs[3];
        index++;

        // ***************************************
        //      Смещение на y - 1
        // ***************************************
        mat->col[index] = realIndex - nx;
        mat->value[index] = coeffs[2];
        index++;

        // ***************************************
        //      Смещение на x - 1
        // ***************************************
        mat->col[index] = realIndex - 1;
        mat->value[index] = coeffs[1];
        index++;

        // Отсутствие смещения
        mat->col[index] = realIndex;
         mat->value[index] = coeffs[0];
        index++;

        // ***************************************
        //      Смещение на x + 1
        // ***************************************
        mat->col[index] = realIndex + 1;
        mat->value[index] = coeffs[1];
        index++;


        // ***************************************
        //      Смещение на y + 1
        // ***************************************
        mat->col[index] = realIndex + nx;
        mat->value[index] = coeffs[2];
        index++;

        // ***************************************
        //      Смещение на z + 1
        // ***************************************
        mat->col[index] = realIndex + nx*ny;
        mat->value[index] = coeffs[3];
        index++;

        k++;
        mat->rowIndex[k] = mat->rowIndex[k-1] + 7;

      }

    }
  }
}

void createImplicitSpMat(SpMatrix* mat, TYPE* coeffs, int nx, int ny, int nz) {
  int index = 0, k = 0;
  int shiftIndexX, shiftIndexY, shiftIndexZ;
  int shift;
  mat->rowIndex[0] = 0;

  int realIndex;
  for (int z = 0; z < nz; z++) {
    if (z == 0)
      shiftIndexZ = nx*ny;
    else if (z == nz - 1)
      shiftIndexZ = -nx*ny;
    else
      shiftIndexZ = 0;
    for (int y = 0; y < ny; y++) {
      if (y == 0)
        shiftIndexY = nx;
      else if (y == ny - 1)
        shiftIndexY = -nx;
      else
        shiftIndexY = 0;

      for (int x = 0; x < nx; x++) {
        if (x == 0)
          shiftIndexX = 1;
        else if (x == nx - 1)
          shiftIndexX = -1;
        else
          shiftIndexX = 0;

        shift = shiftIndexZ + shiftIndexY + shiftIndexX;
        realIndex = x + y*nx + z*nx*ny + shift;

        // ***************************************
        //      Смещение на z - 1
        // ***************************************
        mat->col[index] = realIndex - nx*ny;
        mat->value[index] = coeffs[2];
        index++;

        // ***************************************
        //      Смещение на y - 1
        // ***************************************
        mat->col[index] = realIndex - nx;
        mat->value[index] = coeffs[1];
        index++;

        // ***************************************
        //      Смещение на x - 1
        // ***************************************
        mat->col[index] = realIndex - 1;
        mat->value[index] = coeffs[0];
        index++;

        // ***************************************
        //      Смещение на x + 1
        // ***************************************
        mat->col[index] = realIndex + 1;
        mat->value[index] = coeffs[0];
        index++;


        // ***************************************
        //      Смещение на y + 1
        // ***************************************
        mat->col[index] = realIndex + nx;
        mat->value[index] = coeffs[1];
        index++;

        // ***************************************
        //      Смещение на z + 1
        // ***************************************
        mat->col[index] = realIndex + nx*ny;
        mat->value[index] = coeffs[2];
        index++;

        k++;
        mat->rowIndex[k] = mat->rowIndex[k - 1] + 6;

      }

    }
  }
}