//
// Created by kirill on 19.03.17.
//

#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include <stdbool.h>

#include <sys/utsname.h>
#include <fcntl.h>
#include <unistd.h>

#include "multMV.h"
#include "createSpMat.h"
#include "parser.h"
#include "utils/ts.h"

const char pathSetting[] = "../../../../initial/setting3.ini";
const char pathFunction[] = "../../../../initial/function3.txt";
const char pathResult3D[] = "../../../../result/Kirill/implicit3D_3.txt";

//const char pathSetting[] = "setting3.ini";
//
//const char pathResult3D[] = "res.txt";
//const char pathFunction[] = "function3.txt";

#define EPS 1e-10

#define SHIFT 1
#define RESERVE SHIFT*2

#define IND(x,y,z) ((x) + (y)*NX + (z)*NX*NY)

bool dist(TYPE *x1, TYPE *x2, int N) {
  for (int i = 0; i < N; i++)
    if (fabs(x1[i] - x2[i]) > EPS)
      return false;
  return true;
}

int main() {
  int fd = open("result.txt", O_CREAT | O_RDWR, 0666);
  dup2(fd, 1);


  const size_t len=80;
  char nameHost[len];
  gethostname(nameHost, len);
  printf("name host - %s\n", nameHost);
  fflush(stdout);

  SpMatrix mat;
  Setting setting;
  TYPE coeffs[3];
  TYPE *u;
  int error;
  int NX, NY, NZ;

  error = readSetting(pathSetting, &setting);

  if (error != OK) return error;

  NX = setting.NX + RESERVE;
  NY = setting.NY + RESERVE;
  NZ = setting.NZ + RESERVE;

  int dim = NX*NY*NZ;
  u = (TYPE *)calloc(dim, sizeof(TYPE));

  error = readFunction(pathFunction, u, NX, NY, NZ, SHIFT);

  if (error != OK) return error;

  size_t sizeTime = (size_t) ((setting.TFINISH - setting.TSTART) / setting.dt);

  #if ENABLE_PARALLEL
    printf("PARALLEL VERSION!\n");
    omp_set_num_threads(1);
#endif
  printf("TimeSize -\t%lu\n", sizeTime);

  TYPE dx = fabs(setting.XSTART - setting.XEND) / setting.NX;
  TYPE dy = fabs(setting.YSTART - setting.YEND) / setting.NY;
  TYPE dz = fabs(setting.ZSTART - setting.ZEND) / setting.NZ;

  coeffs[0] = setting.dt*setting.SIGMA/(dx*dx);
  coeffs[1] = setting.dt*setting.SIGMA/(dy*dy);
  coeffs[2] = setting.dt*setting.SIGMA/(dz*dz);

  int nonZero = NX*NY*NZ*6;

  initSpMat(&mat, nonZero, dim);
  createImplicitSpMat(&mat, coeffs, NX, NY, NZ);

//  printSpMat(mat);

  TYPE *x1 = (TYPE *) malloc(sizeof(TYPE) * dim);
  TYPE *x2 = (TYPE *) malloc(sizeof(TYPE) * dim);
  TYPE *tmp;

  int countIter = 0;

  TYPE k = 1.0/(1.0 + 2.0*setting.SIGMA*setting.dt*(1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz)));

  // ОСНОВНЫЕ ВЫЧИСЛЕНИЯ

  double t0 = omp_get_wtime();
  for (int t = 1; t <= sizeTime; t++) {
    memcpy(x1, u, dim*sizeof(TYPE));

    do {
      multMV(x2, mat, x1,0,0,0);


      for (int z = 1; z < NZ-1; z++) {
        for (int y = 1; y < NY-1; y++) {
          for (int x = 1; x < NX-1; x++) {
            x2[IND(x,y,z)] = (u[IND(x,y,z)] + x2[IND(x,y,z)]) * k;
            if (x==1)
              x2[IND(x-1,y,z)]=x2[IND(x,y,z)];
            if (x==NX-2)
              x2[IND(x+1,y,z)]=x2[IND(x,y,z)];
            if (y==1)
              x2[IND(x,y-1,z)]=x2[IND(x,y,z)];
            if (y==NY-2)
              x2[IND(x,y+1,z)]=x2[IND(x,y,z)];
            if (z==1)
              x2[IND(x,y,z-1)]=x2[IND(x,y,z)];
            if (z==NZ-2)
              x2[IND(x,y,z+1)]=x2[IND(x,y,z)];
          }
        }
      }

      tmp = x1;
      x1 = x2;
      x2 = tmp;

      countIter++;

    } while (!dist(x1, x2, dim));

    memcpy(u, x1, dim*sizeof(TYPE));

  }

  //  *******************

  double t1 = omp_get_wtime();
  double diffTime = t1 - t0;
  printf("Count Iteration -\t%d\n", countIter);
  printf("Time -\t%.3lf\n", diffTime);
  writeFunction3D(pathResult3D, u, NX, NY, NZ, SHIFT);

  freeSpMat(&mat);
  free(u);
  free(x1);
  free(x2);
  return 0;
}