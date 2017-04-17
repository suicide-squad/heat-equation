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
#include "sgpu.h"

#define ROOT 0

const char pathSetting[] = "../../../../initial/setting.ini";
const char pathFunction[] = "../../../../initial/function.txt";
const char pathResult1D[] = "../../../../result/Kirill/euler3D_MPI.txt";
const char pathResult3D[] = "../../../../result/Kirill/result_MPI.txt";


int main(int argc, char **argv) {
  int sizeP, rankP;
  size_t sizeTime;

  MPI_Status stats[8];
  double t0, t1;

  int blockYP, blockZP;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankP);
  if (rankP == ROOT) printf("Запуск на %d процессах\n", sizeP);

  SpMatrix mat;
  size_t dim;
  double coeffs[4];
  double* u = NULL, *u_chunk = NULL, *un_chunk;

  size_t NX, NY, NZ, NYr, NZr;

  if (rankP == ROOT) {
    int error;
    Setting setting;
    error = readSetting(pathSetting, &setting);

    if (error != OK) return error;

    dim = (setting.NX + 2) * (setting.NY+2) * (setting.NZ+2);
    u = (double *) malloc(sizeof(double) * dim);
    memset(u, 0, dim * sizeof(double));

    error = readFunction(pathFunction, u, setting.NX+2, setting.NY+2, setting.NZ+2);

    if (error != OK) return error;

    sizeTime = (size_t) ((setting.TFINISH - setting.TSTART) / setting.dt);

    #if ENABLE_PARALLEL
      printf("PARALLEL VERSION 2.0!\n");
    #endif

    printf("TimeSize -\t%lu\n", sizeTime);

    double dx = fabs(setting.XSTART - setting.XEND) / setting.NX;
    double dy = fabs(setting.YSTART - setting.YEND) / setting.NY;
    double dz = fabs(setting.ZSTART - setting.ZEND) / setting.NZ;

    coeffs[0] = 1.0 - 2.0 * setting.dt * setting.SIGMA * (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz));
    coeffs[1] = setting.dt * setting.SIGMA / (dx * dx);
    coeffs[2] = setting.dt * setting.SIGMA / (dy * dy);
    coeffs[3] = setting.dt * setting.SIGMA / (dz * dz);

    NX = setting.NX + 2;
    NY = setting.NY + 2;
    NZ = setting.NZ + 2;
  }

  MPI_Bcast(&sizeTime, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(coeffs, 4, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&NX, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&NY, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&NZ, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);

// TODO:
//  Правильно получить сетку по YZ
//  Ниже представлен костыль на случай 1 и 4 процессов
  if (sizeP != 1) {
    blockYP = (sizeP / 2);
    blockZP = (sizeP / 2);
  } else {
    blockYP = 1;
    blockZP = 1;
  }
  NYr = (NY-2) / blockYP;
  NZr = (NZ-2) / blockZP;

//  СОЗДАНИЕ ДЕКАРДОВОЙ ТОПОЛОГИИ
  const int DIM_CART = 2;
  MPI_Comm gridComm;

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

  u_chunk = (double *)malloc(sizeof(double)*NX*(NYr+2)*(NZr+2));
  un_chunk = (double *)malloc(sizeof(double)*NX*(NYr+2)*(NZr+2));
  double *tmp;

  //  SCATTER
  scatter_by_block(u, u_chunk, NX, NY, NYr, NZr, gridComm);

  int dimChunk = NX*(NYr+2)*(NZr+2);
  int nonZero = dimChunk*7;

  initSpMat(&mat, nonZero, dimChunk);
  createExplicitSpMatV2(&mat, coeffs, NX, NYr + 2, NZr + 2);

  double *bufferRightYSend  = (double *)malloc(sizeof(double)*NX*(NZr+2));
  double *bufferRightYRecv  = (double *)malloc(sizeof(double)*NX*(NZr+2));
  double *bufferLeftYSend   = (double *)malloc(sizeof(double)*NX*(NZr+2));
  double *bufferLeftYRecv   = (double *)malloc(sizeof(double)*NX*(NZr+2));

  double *bufferDownZSend   = (double *)malloc(sizeof(double)*NX*(NYr+2));
  double *bufferDownZRecv   = (double *)malloc(sizeof(double)*NX*(NYr+2));
  double *bufferTopZSend    = (double *)malloc(sizeof(double)*NX*(NYr+2));
  double *bufferTopZRecv    = (double *)malloc(sizeof(double)*NX*(NYr+2));

  int rank_left, rank_right, rank_down, rank_top;
  MPI_Cart_shift(gridComm, 1, 1, &rank_left, &rank_right);
  MPI_Cart_shift(gridComm, 0, 1, &rank_down, &rank_top);

  //printf("rank - %d; left %d; right %d; top %d; down %d\n", rankP, rank_left, rank_right, rank_top, rank_down);

  if (rankP == ROOT) {
    printf("START!\n\n");
    t0 = omp_get_wtime();
  }

  // ОСНОВНЫЕ ВЫЧИСЛЕНИЯ

  for (int t = 1; t <= sizeTime; t++) {
    //  ОБМЕН ГРАНИЦ ПО Y И Z

    //    Передача влево по Y
    pack(bufferLeftYSend, u_chunk, NX, NYr+2, NZr+2, Y_LEFT_SEND);
    MPI_Sendrecv(bufferLeftYSend,  NX*(NZr+2), MPI_DOUBLE, rank_left, 0,
                 bufferRightYRecv, NX*(NZr+2), MPI_DOUBLE, rank_right, 0, gridComm, &stats[0]);
    unpack(bufferRightYRecv, u_chunk, NX, NYr+2, NZr+2, Y_RIGHT_RECV);
    //printf("rank %d, error %d, tag %d, source %d\n",cartrank, stats[0].MPI_ERROR, stats[0].MPI_TAG, stats[0].MPI_SOURCE);

    //    Передача вправо по Y
    pack(bufferRightYSend, u_chunk, NX, NYr+2, NZr+2, Y_RIGHT_SEND);
    MPI_Sendrecv(bufferRightYSend, NX*(NZr+2), MPI_DOUBLE, rank_right, 1,
                 bufferLeftYRecv,  NX*(NZr+2), MPI_DOUBLE, rank_left, 1, gridComm, &stats[1]);
    unpack(bufferLeftYRecv, u_chunk, NX, NYr+2, NZr+2, Y_LEFT_RECV);

    //    Передача вниз по Z
    pack(bufferDownZSend, u_chunk, NX, NYr+2, NZr+2, Z_DOWN_SEND);
    MPI_Sendrecv(bufferDownZSend, NX*(NYr+2), MPI_DOUBLE, rank_down, 2,
                 bufferTopZRecv,  NX*(NYr+2), MPI_DOUBLE, rank_top, 2, gridComm, &stats[2]);
    unpack(bufferTopZRecv, u_chunk, NX, NYr+2, NZr+2, Z_TOP_RECV);

    //    Передача вверх по Z
    pack(bufferTopZSend, u_chunk, NX, NYr+2, NZr+2, Z_TOP_SEND);
    MPI_Sendrecv(bufferTopZSend,  NX*(NYr+2), MPI_DOUBLE, rank_top, 3,
                 bufferDownZRecv, NX*(NYr+2), MPI_DOUBLE, rank_down, 3, gridComm, &stats[3]);
    unpack(bufferDownZRecv, u_chunk, NX, NYr+2, NZr+2, Z_DOWN_RECV);

    //  --------------------------------------

    multMV(&un_chunk, mat, u_chunk);

    tmp = u_chunk;
    u_chunk = un_chunk;
    un_chunk = tmp;

    //  *******************
  }
  if (rankP == ROOT) {
    printf("FINISH!\n\n");
    t1 = omp_get_wtime();
  }

  //        GATHER
  gather_by_block(u, u_chunk, NX, NY, NYr, NZr, gridComm);

  if (rankP == ROOT) {
    double diffTime = t1 - t0;
    printf("Time -\t%.3lf\n", diffTime);
    writeFunction1D(pathResult1D, u, NX);
    writeFunction3D(pathResult3D, u, NX, NY, NZ);

    printf("DONE!!!\n");
    free(u);
  }
  free(bufferRightYSend);
  free(bufferRightYRecv);
  free(bufferLeftYSend);
  free(bufferLeftYRecv);

  free(bufferDownZSend);
  free(bufferDownZRecv);
  free(bufferTopZSend);
  free(bufferTopZRecv);
  free(un_chunk);
  free(u_chunk);
  freeSpMat(&mat);

  MPI_Finalize();
  return 0;
}
