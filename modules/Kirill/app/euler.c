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

    dim = (setting.NX + 2) * setting.NY * setting.NZ;
    u = (double *) malloc(sizeof(double) * dim);
    memset(u, 0, dim * sizeof(double));

    error = readFunction(pathFunction, &u, dim, setting.NX);

    if (error != OK) return error;

    sizeTime = (size_t) ((setting.TFINISH - setting.TSTART) / setting.dt);

    #if ENABLE_PARALLEL
      printf("PARALLEL VERSION 2.0!\n");
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
  NYr = NY / blockYP;
  NZr = NZ / blockZP;

//  СОЗДАНИЕ ДЕКАРДОВОЙ ТОПОЛОГИИ
  const int DIM_CART = 2;
  MPI_Comm gridComm;  	// Коммуникатор в виде квадратной решетки

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

//  TODO:
//  При кубиках на границах достаточно выделять на + 1 (соседей с одной стороны нет)
//  Если ранковый кубик не на границе, то нужно выделять по 2 соседа
// ЛОЖЬ!
//  Всегда нужно знать соседей, при границах зацикливание!!!
  int dim_chunk = NX*(NYr+2)*(NZr+2);

  int nonZero = dim_chunk*7;

  initSpMat(&mat, nonZero, dim_chunk);
  createExplicitSpMat(&mat, coeffs, dim_chunk, NX, NX*(NYr+2));

  double *bufferRightYSend  = (double *)malloc(sizeof(double)*NX*(NZr+2));
  double *bufferRightYRecv  = (double *)malloc(sizeof(double)*NX*(NZr+2));
  double *bufferLeftYSend   = (double *)malloc(sizeof(double)*NX*(NZr+2));
  double *bufferLeftYRecv   = (double *)malloc(sizeof(double)*NX*(NZr+2));

  double *bufferDownZSend   = (double *)malloc(sizeof(double)*NX*(NYr+2));
  double *bufferDownZRecv   = (double *)malloc(sizeof(double)*NX*(NYr+2));
  double *bufferTopZSend    = (double *)malloc(sizeof(double)*NX*(NYr+2));
  double *bufferTopZRecv    = (double *)malloc(sizeof(double)*NX*(NYr+2));

  int rank_recv, rank_send;

  if (rankP == ROOT) {
    printf("START!\n\n");
    t0 = omp_get_wtime();
  }

  // ОСНОВНЫЕ ВЫЧИСЛЕНИЯ

  for (int t = 1; t <= 1; t++) {
    //  ОБМЕН ГРАНИЦ ПО Y И Z

    //    Передача влево по Y
    MPI_Cart_shift(gridComm, 1, -1, &rank_recv, &rank_send);
    pack(bufferLeftYSend, u_chunk, NX, NYr+2, NZr+2, Y_LEFT_SEND);
    MPI_Sendrecv(bufferLeftYSend,  NX*(NZr+2), MPI_DOUBLE, rank_send, 0,
                 bufferRightYRecv, NX*(NZr+2), MPI_DOUBLE, rank_recv, 0, gridComm, &stats[0]);
    unpack(bufferRightYRecv, u_chunk, NX, NYr+2, NZr+2, Y_RIGHT_RECV);

    //    Передача вправо по Y
    MPI_Cart_shift(gridComm, 1, 1, &rank_recv, &rank_send);
    pack(bufferRightYSend, u_chunk, NX, NYr+2, NZr+2, Y_RIGHT_SEND);
    MPI_Sendrecv(bufferRightYSend, NX*(NZr+2), MPI_DOUBLE, rank_send, 1,
                 bufferLeftYRecv,  NX*(NZr+2), MPI_DOUBLE, rank_recv, 1, gridComm, &stats[1]);
    unpack(bufferLeftYRecv, u_chunk, NX, NYr+2, NZr+2, Y_LEFT_RECV);

    //    Передача вниз по Z
    MPI_Cart_shift(gridComm, 0, -1, &rank_recv, &rank_send);
    pack(bufferDownZSend, u_chunk, NX, NYr+2, NZr+2, Z_DOWN_SEND);
    MPI_Sendrecv(bufferDownZSend, NX*(NYr+2), MPI_DOUBLE, rank_send, 2,
                 bufferTopZRecv,  NX*(NYr+2), MPI_DOUBLE, rank_recv, 2, gridComm, &stats[2]);
    unpack(bufferTopZRecv, u_chunk, NX, NYr+2, NZr+2, Z_TOP_RECV);

    //    Передача вверх по Z
    MPI_Cart_shift(gridComm, 0, 1, &rank_recv, &rank_send);
    pack(bufferTopZSend, u_chunk, NX, NYr+2, NZr+2, Z_TOP_SEND);
    MPI_Sendrecv(bufferTopZSend,  NX*(NYr+2), MPI_DOUBLE, rank_send, 3,
                 bufferDownZRecv, NX*(NYr+2), MPI_DOUBLE, rank_recv, 3, gridComm, &stats[3]);
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
    writeFunction1D(pathResult1D, u, NX - 2);
    writeFunction3D(pathResult3D, u, dim, NX - 2);

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
