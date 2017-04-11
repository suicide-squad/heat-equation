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
const char pathResult32D[] = "../../../../result/Kirill/ter.txt";

typedef enum {
  YRIGHT,
  YLEFT,
  ZTOP,
  ZDOWN
} op;

void setex(double *ex, double *u, int NX, int NY, int NZ, op which) {
  switch ( which ) {
//    Y+1
    case YRIGHT: {
      for (int z = 0; z < NZ; z++)
        for (int x = 0; x < NX; x++)
          ex[x + z*NX] = u[x + 0*NX + z*NX*NY];
      break;
    }
//    Y-1
    case YLEFT: {
      for (int z = 0; z < NZ; z++)
        for (int x = 0; x < NX; x++)
          ex[x + z*NX] = u[x + (NY-1)*NX + z*NX*NY];
      break;
    }
      //    Z+1
    case ZTOP: {
      for (int y = 0; y < NY; y++)
        for (int x = 0; x < NX; x++)
          ex[x + y*NX] = u[x + y*NX + (NZ-1)*NX*NY];
      break;
    }
      //    Z-1
    case ZDOWN: {
      for (int y = 0; y < NY; y++)
        for (int x = 0; x < NX; x++)
          ex[x + y*NX] = u[x + y*NX + 0*NX*NY];
      break;
    }
  }
  return;
}

void unpack (double *ex, double *u, int NX, int NY, int NZ, op which) {
  switch ( which ) {
//    Y+1
    case YRIGHT: {
      for (int z = 0; z < NZ; z++)
        for (int x = 0; x < NX; x++)
          u[x + 0*NX + z*NX*NY] = ex[x + z*NX];
      break;
    }
//    Y-1
    case YLEFT: {
      for (int z = 0; z < NZ; z++)
        for (int x = 0; x < NX; x++)
          u[x + (NY-1)*NX + z*NX*NY] = ex[x + z*NX];
      break;
    }
      //    Z+1
    case ZTOP: {
      for (int y = 0; y < NY; y++)
        for (int x = 0; x < NX; x++)
          u[x + y*NX + (NZ-1)*NX*NY] = ex[x + y*NX];
      break;
    }
      //    Z-1
    case ZDOWN: {
      for (int y = 0; y < NY; y++)
        for (int x = 0; x < NX; x++)
          u[x + y*NX + 0*NX*NY] = ex[x + y*NX];
      break;
    }
  }
  return;
}

void scatter_by_block(double *function,
                      double *functionRank,
                      int NX,
                      int NY,
                      int NZ,
                      int NYr,
                      int NZr,
                      MPI_Comm gridComm) {
  int ndim = 0;
  MPI_Cartdim_get(gridComm, &ndim);
  int dims[ndim], periods[ndim], gridCoords[ndim];

  MPI_Cart_get(gridComm, ndim, dims, periods, gridCoords);
  MPI_Comm colComm, rowComm;

  int remain_dims[ndim];
  remain_dims[0] = 1; remain_dims[1] = 0;
  MPI_Cart_sub(gridComm, remain_dims, &colComm);

  int rank;
  MPI_Cart_rank(gridComm, gridCoords,&rank);
  double *functionZ = NULL;
  printf("NZr-%d;NY-%d;NX-%d  r=%d\n",NZr,NY,NX, rank);

  if (gridCoords[1] == 0) {
    functionZ = (double *)malloc(sizeof(double)*NX*NY*(NZr+2));
    MPI_Scatter(function, NX*NY*(NZr), MPI_DOUBLE, functionZ + NX*NY, NX*NY*(NZr), MPI_DOUBLE, 0, colComm);
  }

  printf("scatter rows %d\n", rank);
  remain_dims[0] = 0; remain_dims[1] = 1;
  MPI_Cart_sub(gridComm, remain_dims, &rowComm);

  for (int i = 1; i < NZr + 1; i++)
    MPI_Scatter(functionZ + i*NX*NY, NX*(NYr), MPI_DOUBLE, functionRank + NX + i*NX*(NYr), NX*NYr, MPI_DOUBLE, 0, rowComm);

  if (gridCoords[1] == 0) free(functionZ);
}


void gather_by_block(double *function,
                     double *functionRank,
                     int NX,
                     int NY,
                     int NZ,
                     int NYr,
                     int NZr,
                     MPI_Comm gridComm) {
  int ndim = 0;
  MPI_Cartdim_get(gridComm, &ndim);
  int dims[ndim], periods[ndim], gridCoords[ndim];

  MPI_Cart_get(gridComm, ndim, dims, periods, gridCoords);
  MPI_Comm colComm, rowComm;

  int remain_dims[ndim];
  remain_dims[0] = 1; remain_dims[1] = 0;
  MPI_Cart_sub(gridComm, remain_dims, &colComm);
  remain_dims[0] = 0; remain_dims[1] = 1;
  MPI_Cart_sub(gridComm, remain_dims, &rowComm);

  double *functionZ;
  if (gridCoords[1] == 0) {
    functionZ = (double *) malloc(sizeof(double) * NX * NY *NZr);
  }

  for (int i = 0; i < NZr; i++)
    MPI_Gather(functionRank + NX + (i+1)*(NYr)*NX, NX*(NYr), MPI_DOUBLE, functionZ + i*NX*NY, NX*(NYr), MPI_DOUBLE, 0, rowComm);
  printf("gather cols!\n");


  if (gridCoords[1] == 0) {
    MPI_Gather(functionZ,
               NX * NY * (NZr),
               MPI_DOUBLE,
               function,
               NX * NY * (NZr),
               MPI_DOUBLE,
               0,
               colComm);
    free(functionZ);
  }
}


int main(int argc, char **argv) {
  int sizeP, rankP;
  size_t sizeTime;
  MPI_Status status;
  MPI_Request reqs[8];
  MPI_Status stats[8];
  double t0, t1;

  int blockYP, blockZP;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankP);

  MPI_Comm gridComm;  	// Коммуникатор в виде квадратной решетки

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
  MPI_Bcast(coeffs, 4, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&NX, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&NY, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&NZ, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);

// TODO:
//  Правильно получить сетку по YZ
//  Ниже представлен костыль на случай 4 процессов
  if (sizeP != 1) {
    blockYP = (sizeP / 2);
    blockZP = (sizeP / 2);
  } else {
    blockYP = 1;
    blockZP = 1;
  }
  NYr = NY / blockYP;
  NZr = NZ / blockZP;

  const int DIM_CART = 2;
//  размер каждой размерности
  int dims[DIM_CART];   int periods[DIM_CART];
  int gridCoords[DIM_CART], cartrank;

  dims[0] = blockZP; dims[1] = blockYP;
//  наличие циклов в каждой размерности
  periods[0] = 1; periods[1] = 1;
//  разрешение системе менять номера процессов
  int reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, DIM_CART, dims, periods, reorder, &gridComm);

  // Координаты текущего процесса в процессной решетке и ранг в решётке
  MPI_Comm_rank(gridComm, &cartrank);
  // Определение координат процесса в решетке
  MPI_Cart_coords(gridComm, rankP, DIM_CART, gridCoords);

  functionRank = (double *)malloc(sizeof(double)*NX*(NYr+2)*(NZr+2));

  scatter_by_block(function, functionRank, NX, NY, NZ, NYr, NZr, gridComm);

  printf("scatter finish!\n");

//  TODO:
//  При кубиках на границах достаточно выделять на + 1 (соседей с одной стороны нет)
//  Если ранковый кубик не на границе, то нужно выделять по 2 соседа
// ЛОЖЬ!
//  Всегда нужно знать соседей, при границах зацикливание!!!
  size_t dimPart = NX*(NYr)*(NZr);
  //  SCATTER

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

  double *bufferRightY = (double *)malloc(sizeof(double)*NX*NZr);
  double *bufferLeftY = (double *)malloc(sizeof(double)*NX*NZr);
  double *bufferTopZ = (double *)malloc(sizeof(double)*NX*NYr);
  double *bufferDownZ = (double *)malloc(sizeof(double)*NX*NYr);

  int nextY, prevY, nextZ, prevZ;

  if (rankP == ROOT) t0 = omp_get_wtime();

  // ОСНОВНЫЕ ВЫЧИСЛЕНИЯ

//  for (int t = 1; t <= 1; t++) {
//////    TODO:
//////    ТУТ должен быть обмен границ!
////    Отправка граничных условий
//
//
////  принять и распаковать границы по Z
//    unpack(bufferDownZ, functionRank, NX, NYr, NZr, ZDOWN);
//    unpack(bufferTopZ, functionRank, NX, NYr, NZr, ZTOP);
//    multMV(&nextFunction, mat, functionRank);
//
//    tmp = functionRank;
//    function = nextFunction;
//    nextFunction = tmp;
//
//    //  *******************
//  }
  if (rankP == ROOT) t1 = omp_get_wtime();

  printf("rank = %d\n",rankP);
//        GATHER
//  Нулевой процесс собирает сразу, остальные отсылают на сборку

  gather_by_block(function, functionRank, NX, NY, NZ, NYr, NZr, gridComm);

  if (rankP == ROOT) {
    double diffTime = t1 - t0;
    printf("Time -\t%.3lf\n", diffTime);
    writeFunction1D(pathResult1D, function, NX - 2);
    writeFunction3D(pathResult3D, function, dim, NX - 2);

    free(function);
  }
  free(bufferRightY);
  free(bufferLeftY);
  free(bufferDownZ);
  free(bufferTopZ);
  freeSpMat(&mat);
  free(nextFunction);
  free(functionRank);

  MPI_Finalize();
  return 0;
}
