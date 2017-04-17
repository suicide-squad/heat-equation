//
// Created by kirill on 19.03.17.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>

#include "sp_mat.h"
#include "parser.h"

const char pathSetting[] = "../../../../initial/setting.ini";
const char pathFunction[] = "../../../../initial/function.txt";
const char pathResult[] = "../../../../result/Kirill/runge3D.txt";
const char pathResult3D[] = "../../../../result/Kirill/result_runge3D.txt";

int main() {
  SpMatrix A, B, C;
  Setting setting;
  double coeffs[4];
  double* u;
  int error;
  error = readSetting(pathSetting, &setting);

  if (error != OK) return error;

  int NX = (setting.NX + 2);
  int NY = (setting.NY + 2);
  int NZ = (setting.NZ + 2);

  size_t dim = NX*NY*NZ;
  u = (double *)malloc(sizeof(double)*dim);
  memset(u, 0, dim* sizeof(double));

  error = readFunction(pathFunction, u, NX, NY, NZ);

  if (error != OK) return error;

  size_t sizeTime = (size_t)((setting.TFINISH - setting.TSTART) / setting.dt);


  #if ENABLE_PARALLEL
  printf("PARALLEL VERSION!\n");
  #endif
  printf("TimeSize -\t%lu\n", sizeTime);

  double hX = fabs(setting.XSTART - setting.XEND) / setting.NX;
  double hY = fabs(setting.YSTART - setting.YEND) / setting.NY;
  double hZ = fabs(setting.ZSTART - setting.ZEND) / setting.NZ;

  coeffs[0] = setting.SIGMA/(hX*hX);
  coeffs[2] = setting.SIGMA/(hY*hY);
  coeffs[3] = setting.SIGMA/(hZ*hZ);
  coeffs[1] = - 2.0*setting.SIGMA*(coeffs[0] + coeffs[2] + coeffs[3]);

  size_t nonZero = dim*7;

  initSpMat(&A, nonZero, dim);
  createExplicitSpMatV2(&A, coeffs, NX, NY, NZ);

  coeffs[0] = setting.dt*coeffs[0]*0.5;
  coeffs[2] = setting.dt*coeffs[2]*0.5;
  coeffs[3] = setting.dt*coeffs[3]*0.5;
  coeffs[1] = 1.0 - 2.0*(coeffs[0] + coeffs[2] + coeffs[3]);

  initSpMat(&B, nonZero, dim);
  createExplicitSpMatV2(&B, coeffs, NX, NY, NZ);

  coeffs[0] = coeffs[0]*2.0;
  coeffs[2] = coeffs[2]*2.0;
  coeffs[3] = coeffs[3]*2.0;
  coeffs[1] = 1.0 - 2.0*(coeffs[0] + coeffs[2] + coeffs[3]);

  initSpMat(&C, nonZero, dim);
  createExplicitSpMatV2(&C, coeffs, NX, NY, NZ);

  //  printSpMat(A);

  double *un = (double *)malloc(sizeof(double)*dim);
  double* k1 = (double*)malloc(sizeof(double)*dim);
  double* k2 = (double*)malloc(sizeof(double)*dim);
  double* k3 = (double*)malloc(sizeof(double)*dim);
  double* k4 = (double*)malloc(sizeof(double)*dim);
  double h = setting.dt/6.0;
  double *tmp;

  double t0 = omp_get_wtime();

  // ОСНОВНЫЕ ВЫЧИСЛЕНИЯ

  for (int t = 1; t <= 1; t++) {
//    ОБНОВИТЬ ПЕРЕДАЧУ ГРАНИЦ!!!

    // k1 = A*U
    multMV(&k1, A, u);

    // k2 = B*k1
    multMV(&k2, B, k1);

    // k3 = B*k2
    multMV(&k3, B, k2);

    // k4 = C*k3
    multMV(&k4, C, k3);

    // UNext = U + (k1 + k2*2 + k3*2 + k4)*h;
    sumV(&un, u, k1, k2, k3, k4, dim, h);

    tmp = u;
    u = un;
    un = tmp;
  }

  //  *******************

  double t1 = omp_get_wtime();
  double diffTime = t1 - t0;
  printf("Time -\t%.3lf\n", diffTime);
  writeFunction1D(pathResult, u, NX);
  writeFunction3D(pathResult3D, u, NX, NY, NZ);

  freeSpMat(&A);

  free(u);
  free(un);
  free(k1);
  free(k2);
  free(k3);
  free(k4);
  return 0;
}
