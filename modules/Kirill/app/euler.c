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
const char pathResult[] = "../../../../result/Kirill/euler.txt";

int main() {
  SpMatrix mat;
  Setting setting;
  double coeffs[4];
  double* function;
  int error;
  error = readSetting(pathSetting, &setting);

  if (error != OK) return error;

  size_t dim = (setting.NX + 2)*setting.NY*setting.NZ;
  function = (double *)malloc(sizeof(double)*dim);
  memset(function, 0, dim* sizeof(double));

  error = readFunction(pathFunction, &function, dim, setting.NX);

  if (error != OK) return error;

  size_t sizeTime = (size_t)((setting.TFINISH - setting.TSTART) / setting.dt);

  #if ENABLE_PARALLEL
  printf("PARALLEL VERSION!\n");
  #endif
  printf("TimeSize -\t%lu\n", sizeTime);

  double hX = fabs(setting.XSTART - setting.XEND) / setting.NX;
  double hY = fabs(setting.YSTART - setting.YEND) / setting.NY;
  double hZ = fabs(setting.ZSTART - setting.ZEND) / setting.NZ;

  coeffs[0] = setting.dt*setting.SIGMA/(hX*hX);
  coeffs[1] = 1.0 - 2.0*setting.dt*setting.SIGMA*(1.0/(hX*hX) + 1.0/(hY*hY) + 1.0/(hZ*hZ));
  coeffs[2] = setting.dt*setting.SIGMA/(hY*hY);
  coeffs[3] = setting.dt*setting.SIGMA/(hZ*hZ);

  size_t NZ = setting.NX*setting.NY*setting.NZ*7 + (dim - setting.NX*setting.NY*setting.NZ)*3;
  int NX = (int)setting.NX + 2;
  int NXY = (int)setting.NY*NX;

  initSpMat(&mat, NZ, dim);
  createExplicitSpMat(&mat, coeffs, (int)dim, NX, NXY);

//  printSpMat(mat);

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
  free(function);
  free(nextFunction);
  return 0;
}
