//
// Created by kirill on 7.03.17.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <parser.h>
#include <memory.h>

#include "sp_mat.h"
#include "parser.h"

const char pathSetting[] = "../../../../initial/setting.ini";
const char pathFunction[] = "../../../../initial/function.txt";

int main() {

  Setting setting;
  double* function;
  int error;
  double** MATRIX;
  error = readSetting(pathSetting, &setting);

  if (error != OK) return error;

  size_t dim = (setting.NX+2)*setting.NY*setting.NZ;
  function = (double *)malloc(sizeof(double)*dim);
  memset(function, 0, dim);

  error = readFunction(pathFunction, &function, dim, setting.NX);

  if (error != OK) return error;

  size_t sizeTime = (size_t)((setting.TFINISH - setting.TSTART) / setting.dt);

  printf("%lu\n", sizeTime);

  double hX = fabs(setting.XSTART - setting.XEND) / setting.NX;
  double hY = fabs(setting.YSTART - setting.YEND) / setting.NY;
  double hZ = fabs(setting.ZSTART - setting.ZEND) / setting.NZ;

  double IJK = 1.0 - 2.0*setting.dt*setting.SIGMA*(1.0/(hX*hX) + 1.0/(hY*hY) + 1.0/(hZ*hZ));
  double Ijk = setting.dt*setting.SIGMA/(hX*hX);
  double iJk = setting.dt*setting.SIGMA/(hY*hY);
  double ijK = setting.dt*setting.SIGMA/(hZ*hZ);

  MATRIX = (double **)malloc(sizeof(double*)*dim);
  for (int i = 0; i < dim; i++) {
    MATRIX[i] = (double *) malloc(sizeof(double) * dim);
    memset(MATRIX[i], 0, dim);
  }

  int j,k;
  int NX = (int)setting.NX + 2;
  for (int i = 0; i < dim; i++) {
//    Условие, чтобы учесть граничные условия по X
    if ( i%NX!= 0 && i%NX != NX - 1 ) {
      MATRIX[i][i] = 1.0;
      MATRIX[i][i - 1] = 2.0;
      MATRIX[i][i + 1] = 2.0;
      j = i - i*NX;
//      Цикличность условия по Y
      if (j <= 0)
        MATRIX[i][(dim - j%dim)] = 3.0;
      else
        MATRIX[i][j%dim] = 3.0;
      MATRIX[i][(i + i*NX) % dim] = 3.0;
//    Цикличность условия по Z
      k = j - i*NX*(int)setting.NY;
      if (k <= 0)
        MATRIX[i][dim - k%dim] = 4.0;
      else
        MATRIX[i][k%dim] = 4.0;
      MATRIX[i][(i + i*NX + i*NX*setting.NY) % dim] = 4.0;
    } else
      MATRIX[i][i] = 5.0;

  }

  double *nextFunction = (double *)malloc(sizeof(double)*dim);
  double *tmp;

  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      printf("%.0lf ", MATRIX[i][j]);

    }
    printf("\n");
  }


//  for (int t = 1; t <= sizeTime; t++) {
//    mult(&nextFunction, MATRIX, function, dim);
//
//    tmp = function;
//    function = nextFunction;
//    nextFunction = tmp;
//  }
//
//  writeFunctionX("test.txt", function, dim, (size_t )NX);

  return 0;

}
