//
// Created by kirill on 7.03.17.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>

#include "sp_mat.h"
#include "parser.h"

const char pathSetting[] = "../../../../initial/setting.ini";
const char pathFunction[] = "../../../../initial/function.txt";
const char pathResult[] = "../../../../result/Kirill/result.txt";

int main() {
  SpMatrix mat;
  Setting setting;
  double* function;
  int error;
  error = readSetting(pathSetting, &setting);

  if (error != OK) return error;

  size_t dim = (setting.NX+2)*setting.NY*setting.NZ;
  function = (double *)malloc(sizeof(double)*dim);
  memset(function, 0, dim* sizeof(double));

  error = readFunction(pathFunction, &function, dim, setting.NX);

//  for (int i = 0; i < dim; i++) {
//    printf("%.0lf\n", function[i]);
//  }

  if (error != OK) return error;

  size_t sizeTime = (size_t)((setting.TFINISH - setting.TSTART) / setting.dt);

  #if ENABLE_PARALLEL
  printf("ПАРАЛЛЕЛЬНАЯ ВЕРСИЯ!\n");
  #endif
  printf("Количество шагов по времени -\t%lu\n", sizeTime);

  double hX = fabs(setting.XSTART - setting.XEND) / setting.NX;
  double hY = fabs(setting.YSTART - setting.YEND) / setting.NY;
  double hZ = fabs(setting.ZSTART - setting.ZEND) / setting.NZ;

  double IJK = 1.0 - 2.0*setting.dt*setting.SIGMA*(1.0/(hX*hX) + 1.0/(hY*hY) + 1.0/(hZ*hZ));
  double Ijk = setting.dt*setting.SIGMA/(hX*hX);
  double iJk = setting.dt*setting.SIGMA/(hY*hY);
  double ijK = setting.dt*setting.SIGMA/(hZ*hZ);

  size_t NZ = setting.NX*setting.NY*setting.NZ*7 + (dim - setting.NX*setting.NY*setting.NZ);
  initSpMat(&mat, NZ, dim);
  int j,k;
  int NX = (int)setting.NX + 2;
  int NXY = (int)setting.NY*NX;
  int index = 0;
  mat.rowIndex[0] = 0;

  for (int i = 0; i < dim; i++) {
    if ( i%NX!= 0 && i%NX != NX - 1 ) {
      // Смещение на x + 1 и x - 1 с учётом граничных условий
      // ***************************************
      mat.col[index] = i - 1;
      mat.value[index] = Ijk;
      index++;

      mat.col[index] = i;
      mat.value[index] = IJK;
      index++;

      mat.col[index] = i + 1;
      mat.value[index] = Ijk;
      index++;
      // ***************************************

      // Смещение на y + 1 и y - 1 с учётом цикличности условия
      // ***************************************
      j = i - NX;
      if (j <= 0) {
        mat.col[index] = (int)dim + j;
        mat.value[index] = iJk;
        index++;
      }
      else {
        mat.col[index] = j;
        mat.value[index] = iJk;
        index++;
      }
      mat.col[index] = (i + NX) % NXY;
      mat.value[index] = iJk;
      index++;
      // ***************************************

      // Смещение на z + 1 и z - 1 с учётом цикличности условия
      // ***************************************
      k = i - NXY;
      if (k <= 0) {
        mat.col[index] = (int)dim + k;
        mat.value[index] = ijK;
        index++;
      }
      else {
        mat.col[index] = k;
        mat.value[index] = ijK;
        index++;

      }
      mat.col[index] = (i + NXY)%(int)dim;
      mat.value[index] = ijK;
      index++;
      // ***************************************

      mat.rowIndex[i + 1] = mat.rowIndex[i] + 7;
    } else {
      mat.col[index] = i;
      mat.value[index] = 1.0;
      index++;

      mat.rowIndex[i + 1] = mat.rowIndex[i] + 1;
    }
  }

  double *nextFunction = (double *)malloc(sizeof(double)*dim);
  double *tmp;

  double t0 = omp_get_wtime();

  // ОСНОВНЫЕ ВЫЧИСЛЕНИЯ

  for (int t = 1; t <= sizeTime; t++) {
    multMV(&nextFunction, mat, function);

    tmp = function;
    function = nextFunction;
    nextFunction = tmp;
  }

  //  *******************

  double t1 = omp_get_wtime();
  double diffTime = t1 - t0;
  printf("Time -\t%.3lf\n", diffTime);
  writeFunctionX(pathResult, function, setting.NX);

  freeSpMat(&mat);
  return 0;
}
