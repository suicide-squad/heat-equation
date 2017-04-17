//
// Created by kirill on 19.03.17.
//

#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include <stdbool.h>

#include "sp_mat.h"
#include "parser.h"

const double EPS = 1e-10;

const char pathSetting[] = "../../../../initial/setting.ini";
const char pathFunction[] = "../../../../initial/function.txt";
const char pathResult[] = "../../../../result/Kirill/implicit3D.txt";

bool dist(double *x1, double *x2, size_t N) {
  for (int i = 0; i < N; i++)
    if (fabs(x1[i] - x2[i]) > EPS)
      return false;
  return true;
}

int main() {
  SpMatrix mat;
  Setting setting;
  double coeffs[3];
  double *function;
  int error;
  error = readSetting(pathSetting, &setting);

  if (error != OK) return error;

  size_t dim = (setting.NX + 2) * setting.NY * setting.NZ;
  function = (double *) malloc(sizeof(double) * dim);
  memset(function, 0, dim * sizeof(double));

  error = readFunction(pathFunction, &function, setting.NX, 0, 0);

  if (error != OK) return error;

  size_t sizeTime = (size_t) ((setting.TFINISH - setting.TSTART) / setting.dt);

#if ENABLE_PARALLEL
  printf("PARALLEL VERSION!\n");
#endif
  printf("TimeSize -\t%lu\n", sizeTime);

  double hX = fabs(setting.XSTART - setting.XEND) / setting.NX;
  double hY = fabs(setting.YSTART - setting.YEND) / setting.NY;
  double hZ = fabs(setting.ZSTART - setting.ZEND) / setting.NZ;

  coeffs[0] = setting.dt * setting.SIGMA / (hX * hX);
  coeffs[1] = setting.dt * setting.SIGMA / (hY * hY);
  coeffs[2] = setting.dt * setting.SIGMA / (hZ * hZ);

  size_t NZ = setting.NX*setting.NY*setting.NZ*6;
  int NX = (int) setting.NX + 2;
  int NXY = (int) setting.NY * NX;

  initSpMat(&mat, NZ, dim);
  createImplicitSpMat(&mat, coeffs, (int) dim, NX, NXY);

//  printSpMat(mat);

  double *X1 = (double *) malloc(sizeof(double) * dim);
  double *X2 = (double *) malloc(sizeof(double) * dim);
  double *tmp;

  double t0 = omp_get_wtime();

  // ОСНОВНЫЕ ВЫЧИСЛЕНИЯ

  int countIter = 0;

  double k = 1.0/(1.0 + 2.0*setting.SIGMA*setting.dt*(1.0/(hX*hX) + 1.0/(hY*hY) + 1.0/(hZ*hZ)));

  for (int t = 1; t <= sizeTime; t++) {
    memcpy(X1, function, dim*sizeof(double));

    do {
      multMV(&X2, mat, X1);

      for (int i = 0; i < dim; i++) {
        if (i % NX != 0 && i % NX != NX - 1) {
          X2[i] = (function[i] + X2[i]) * k;
        }
        if (i % NX == 1) {
          X2[i - 1] = X2[i];
        }
        if (i % NX == NX - 1) {
          X2[i] = X2[i - 1];
        }
      }

      tmp = X1;
      X1 = X2;
      X2 = tmp;

      countIter++;

    } while (!dist(X1, X2, dim));

    memcpy(function, X1, dim*sizeof(double));

  }

  //  *******************

  double t1 = omp_get_wtime();
  double diffTime = t1 - t0;
  printf("Time -\t%.3lf\n", diffTime);
  writeFunction1D(pathResult, function, setting.NX);

  freeSpMat(&mat);
  free(function);
  free(X1);
  free(X2);
  return 0;
}