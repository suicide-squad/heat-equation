//
// Created by kirill on 7.03.17.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <mpi.h>
#include <parser.h>

#include "sp_mat.h"
#include "parser.h"

#define ROOT 0

const char pathSetting[] = "../../../../initial/setting.ini";
const char pathFunction[] = "../../../../initial/function.txt";
const char pathResult1D[] = "../../../../result/Kirill/euler3D_MPI.txt";
const char pathResult3D[] = "../../../../result/Kirill/result_MPI.txt";

int main(int argc, char **argv) {
  int sizeP, rankP;
  size_t sizeTime;
  MPI_Status status;
  double t0, t1;


  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankP);

  SpMatrix mat;
  size_t dim;
  double coeffs[4];
  double* function;
  double* functionRank;


  size_t NX, NY, NZ, NYr, NZr;

  if (rankP == ROOT) {
    int error;
    Setting setting;
    error = readSetting(pathSetting, &setting);

    if (error != OK) return error;

    dim = (setting.NX + 2) * setting.NY * setting.NZ;
    function = (double *) malloc(sizeof(double) * dim);
    memset(function, 0, dim * sizeof(double));

    error = readFunction(pathFunction, &function, dim, setting.NX);

    if (error != OK) return error;

    sizeTime = (size_t) ((setting.TFINISH - setting.TSTART) / setting.dt);

    #if ENABLE_PARALLEL
      printf("PARALLEL VERSION!\n");
    #endif

    printf("TimeSize -\t%lu\n", sizeTime);

    double hX = fabs(setting.XSTART - setting.XEND) / setting.NX;
    double hY = fabs(setting.YSTART - setting.YEND) / setting.NY;
    double hZ = fabs(setting.ZSTART - setting.ZEND) / setting.NZ;

    coeffs[0] = setting.dt * setting.SIGMA / (hX * hX);
    coeffs[1] = 1.0 - 2.0 * setting.dt * setting.SIGMA * (1.0 / (hX * hX) + 1.0 / (hY * hY) + 1.0 / (hZ * hZ));
    coeffs[2] = setting.dt * setting.SIGMA / (hY * hY);
    coeffs[3] = setting.dt * setting.SIGMA / (hZ * hZ);

    NX = setting.NX + 2;
    NY = setting.NY;
    NZ = setting.NZ;
  }

  MPI_Bcast(&sizeTime, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&dim, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&NX, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&NY, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&NZ, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);

// TODO:
//  Правильно получить сетку по YZ
//  Ниже представлен костыль на случай 4 процессов
  if (sizeP != 1) {
    NYr = NY / (sizeP / 2);
    NZr = NZ / (sizeP / 2);
  }

//  TODO:
//  При кубиках на границах достаточно выделять на + 1 (соседей с одной стороны нет)
//  Если ранковый кубик не на границе, то нужно выделять по 2 соседа
// ЛОЖЬ!
//  Всегда нужно знать соседей, при границах зацикливание!!!
  size_t dimPart = NX*(NYr + 2)*(NZr + 2);
  functionRank = (double *)malloc(sizeof(double)*dimPart);

  //  SCATTER
  if (rankP == ROOT) {
    for (int rank = 1; rank < sizeP; rank++) {
      for (int z = 1; z < NZr + 1; z++)
        for (int y = 1; y < NYr + 1; y++)
          for (int x = 0; x < NX; x++)
            functionRank[x + y*NX + z*NX*NYr] = function[x + (y + (rank%2)*NYr)*NX + (z + (rank/2)*NZr)*NX*NY];

      MPI_Rsend(functionRank, dimPart, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD);
    }
    for (int z = 1; z < NZr + 1; z++)
      for (int y = 1; y < NYr + 1; y++)
        for (int x = 0; x < NX; x++)
          functionRank[x + y*NX + z*NX*NYr] = function[x + y*NX + z*NX*NY];

  }

  if (rankP != ROOT) {
    MPI_Recv(functionRank, dimPart, MPI_DOUBLE, ROOT, 0, MPI_COMM_WORLD, &status);
  }

//  ДЛЯ ОТЛАДКИ
//  setting.NX = 3;
//  setting.NY = 3;
//  setting.NZ = 3;
//  dim = 45;
//  coeffs[0] = 1;
//  coeffs[1] = 2;
//  coeffs[2] = 3;
//  coeffs[3] = 4;

  size_t NonZero = dimPart*7;
  int NXY = NYr*NX;

  initSpMat(&mat, NonZero, dimPart);
  createExplicitSpMat(&mat, coeffs, dimPart, NX, NXY);

//  printSpMat(mat);

  double *nextFunction = (double *)malloc(sizeof(double)*dimPart);
  double *tmp;

  if (rankP == ROOT) t0 = omp_get_wtime();

  // ОСНОВНЫЕ ВЫЧИСЛЕНИЯ

//  for (int t = 1; t <= sizeTime; t++) {
////    TODO:
////    ТУТ должен быть обмен границ!
//
//    multMV(&nextFunction, mat, functionRank);
//
//    tmp = functionRank;
//    function = nextFunction;
//    nextFunction = tmp;
//  }

  //  *******************

  if (rankP == ROOT) t1 = omp_get_wtime();

//        GATHER
//  Нулевой процесс собирает сразу, остальные отсылают на сборку
  if (rankP == ROOT) {
    for (int z = 1; z < NZr; z++)
      for (int y = 1; y < NYr; y++)
        for (int x = 0; x < NX; x++)
          function[x + (y-1)*NX + (z-1)*NX*NY] = functionRank[x + y*NX + z*NX*NYr];
  }
  else
    MPI_Rsend(functionRank, dimPart, MPI_DOUBLE, ROOT, 0, MPI_COMM_WORLD);
//  Сбор в один вектор от всех процессов
  if (rankP == ROOT) {
    for (int rank = 1; rank < sizeP; rank++) {
      MPI_Recv(functionRank, dimPart, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, &status);
      for (int z = 1; z < NZr; z++)
        for (int y = 1; y < NYr; y++)
          for (int x = 0; x < NX; x++)
            function[x + ((y-1)+(rank%2)*(NYr-2))*NX + ((z-1)+(rank/2)*(NZr-2))*NX*NY] =
                functionRank[x + y*NX + z*NX*NYr];
    }
  }

  if (rankP == ROOT) {
    double diffTime = t1 - t0;
    printf("Time -\t%.3lf\n", diffTime);
    writeFunction1D(pathResult1D, function, NX - 2);
    writeFunction3D(pathResult3D, function, dim, NX - 2);

    free(function);
  }
  freeSpMat(&mat);
  free(nextFunction);
  free(functionRank);

  MPI_Finalize();
  return 0;
}
