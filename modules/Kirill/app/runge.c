//
// Created by kirill on 19.03.17.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

#include <sys/utsname.h>

#include "multMV.h"
#include "parser.h"
#include "sgpu.h"
#include "createSpMat.h"

#include "utils/ts.h"

#define DIM_CART 2
#define SHIFT 4
#define RESERVE SHIFT*2

#define IND(x,y,z) ((x) + (y)*NX + (z)*NX*(NYr + RESERVE))

const char pathSetting[] = "../../../../initial/setting3.ini";
const char pathFunction[] = "../../../../initial/function3.txt";
const char pathResult3D[] = "../../../../result/Kirill/runge3D_3.txt";


int main(int argc, char **argv) {
  int sizeP, rankP;
  size_t sizeTime;
  double t1 = 0.0, t0 = 0.0;
  SpMatrix A, B, C;
  Setting setting;
  TYPE coeffs[4];
  TYPE dt;
  MPI_Status status[4];
  int blockYP = 0, blockZP = 0;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankP);

  const size_t len=80;
  char nameHost[len];
  gethostname(nameHost, len);
  printf("rank - %d name host - %s\n",rankP, nameHost);

  TYPE * u = NULL, *u_chunk = NULL, *un_chunk = NULL;
  int NX, NY, NZ, NYr, NZr;

#if ENABLE_PARALLEL
        omp_set_num_threads(atoi(argv[1]));
        
    if (rankP == ROOT) printf("PARALLEL VERSION! Number of threads - %u\n", omp_get_max_threads());
#endif

  if (rankP == ROOT) {

    int error = readSetting(pathSetting, &setting);
    if (error != OK) return error;

    NX = (setting.NX + RESERVE);
    NY = (setting.NY + RESERVE);
    NZ = (setting.NZ + RESERVE);

    int dim = NX*NY*NZ;
    u = (TYPE *) calloc(dim, sizeof(TYPE));
    error = readFunction(pathFunction, u, NX, NY, NZ, SHIFT);

    if (error != OK) return error;

    sizeTime = (size_t) ((setting.TFINISH - setting.TSTART)/setting.dt);

    printf("TimeSize -\t%lu\n", sizeTime);

    TYPE dx = ABS(setting.XSTART - setting.XEND)/setting.NX;
    TYPE dy = ABS(setting.YSTART - setting.YEND)/setting.NY;
    TYPE dz = ABS(setting.ZSTART - setting.ZEND)/setting.NZ;

    dt = setting.dt;

    coeffs[1] = setting.SIGMA/(dx*dx);
    coeffs[2] = setting.SIGMA/(dy*dy);
    coeffs[3] = setting.SIGMA/(dz*dz);
    coeffs[0] = -2.0*(coeffs[1] + coeffs[2] + coeffs[3]);
  }

  MPI_Bcast(&sizeTime, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(coeffs, 4, MPI_TYPE, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&dt, 1, MPI_TYPE, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&NX, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&NY, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&NZ, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

  //  Определения числа процессов в каждом измерении
  get_blocks(&blockYP, &blockZP, sizeP);

  if (rankP == ROOT) printf("blockY %d blockZ %d\n", blockYP, blockZP);

  NYr = (NY - RESERVE)/blockYP;
  NZr = (NZ - RESERVE)/blockZP;

  MPI_Comm gridComm;

  int dims[DIM_CART], periods[DIM_CART];
  int gridCoords[DIM_CART];

  //  размер каждой размерности
  dims[0] = blockZP; dims[1] = blockYP;
  //  наличие циклов в каждой размерности
  periods[0] = 0; periods[1] = 0;
  //  разрешение системе менять номера процессов
  int reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, DIM_CART, dims, periods, reorder, &gridComm);

  // Определение координат процесса в решетке
  MPI_Cart_coords(gridComm, rankP, DIM_CART, gridCoords);

  int dimChunk = NX*(NYr + RESERVE)*(NZr + RESERVE);
  u_chunk = (TYPE *)malloc(sizeof(TYPE)*dimChunk);
  un_chunk = (TYPE *)malloc(sizeof(TYPE)*dimChunk);

  //  SCATTER
  scatter_by_block(u, u_chunk, NX, NY, NYr, NZr, gridComm, RESERVE);
  // if (rankP == ROOT) printf("Scatter!\n");

  int nonZero = dimChunk*7;
    
    // Определения соседних ранков в декардовой решётке
    int rank_left, rank_right, rank_down, rank_top;
    MPI_Cart_shift(gridComm, 1, -1, &rank_left, &rank_right);
    MPI_Cart_shift(gridComm, 0, 1, &rank_down, &rank_top);
    
    
  initSpMat(&A, nonZero, dimChunk);
  createExplicitSpMatV2R(&A, coeffs, NX, NYr + RESERVE, NZr + RESERVE, rank_left, rank_right, rank_down, rank_top);

  coeffs[1] = dt*coeffs[1]*0.5;
  coeffs[2] = dt*coeffs[2]*0.5;
  coeffs[3] = dt*coeffs[3]*0.5;
  coeffs[0] = 1.0 - 2.0*(coeffs[1] + coeffs[2] + coeffs[3]);

  initSpMat(&B, nonZero, dimChunk);
  createExplicitSpMatV2R(&B, coeffs, NX, NYr + RESERVE, NZr + RESERVE, rank_left, rank_right, rank_down, rank_top);

  coeffs[1] = coeffs[1]*2.0;
  coeffs[2] = coeffs[2]*2.0;
  coeffs[3] = coeffs[3]*2.0;
  coeffs[0] = 1.0 - 2.0*(coeffs[1] + coeffs[2] + coeffs[3]);

  initSpMat(&C, nonZero, dimChunk);
  createExplicitSpMatV2R(&C, coeffs, NX, NYr + RESERVE, NZr + RESERVE, rank_left, rank_right, rank_down, rank_top);

  //  printSpMat(A);

  TYPE* k1 = (TYPE*)malloc(sizeof(TYPE)*dimChunk);
  TYPE* k2 = (TYPE*)malloc(sizeof(TYPE)*dimChunk);
  TYPE* k3 = (TYPE*)malloc(sizeof(TYPE)*dimChunk);
  TYPE* k4 = (TYPE*)malloc(sizeof(TYPE)*dimChunk);
  TYPE h = dt/6.0;
  TYPE *tmp;

  for (int z = 0; z < NZr + RESERVE; z++) {
    for (int y = 0; y < NYr + RESERVE; y++) {
      for (int x = 0; x < NX; x++) {
        if (x == SHIFT - 1) {
          u_chunk[IND(x, y, z)] = u_chunk[IND(x + 1, y, z)];
          u_chunk[IND(x-1, y, z)] = u_chunk[IND(x, y, z)];
          u_chunk[IND(x-2, y, z)] = u_chunk[IND(x-1, y, z)];
          u_chunk[IND(x-3, y, z)] = u_chunk[IND(x-2, y, z)];
        }
        if (x == NX - SHIFT ) {
          u_chunk[IND(x, y, z)] = u_chunk[IND(x -1, y, z)];
          u_chunk[IND(x+1, y, z)] = u_chunk[IND(x, y, z)];
          u_chunk[IND(x+2, y, z)] = u_chunk[IND(x+1, y, z)];
          u_chunk[IND(x+3, y, z)] = u_chunk[IND(x+2, y, z)];
        }
        if (y == SHIFT  - 1) {
          u_chunk[IND(x, y, z)] = u_chunk[IND(x, y + 1, z)];
          u_chunk[IND(x, y-1, z)] = u_chunk[IND(x, y, z)];
          u_chunk[IND(x, y-2, z)] = u_chunk[IND(x, y - 1, z)];
          u_chunk[IND(x, y-3, z)] = u_chunk[IND(x, y - 2, z)];
        }
        if (y == NYr + SHIFT ) {
          u_chunk[IND(x, y, z)] = u_chunk[IND(x, y - 1, z)];
          u_chunk[IND(x, y + 1, z)] = u_chunk[IND(x, y, z)];
          u_chunk[IND(x, y + 2, z)] = u_chunk[IND(x, y + 1, z)];
          u_chunk[IND(x, y + 3, z)] = u_chunk[IND(x, y + 2, z)];
        }
        if (z == SHIFT - 1) {
          u_chunk[IND(x, y, z)] = u_chunk[IND(x, y, z + 1)];
          u_chunk[IND(x, y, z - 1)] = u_chunk[IND(x, y, z)];
          u_chunk[IND(x, y, z - 2)] = u_chunk[IND(x, y, z - 1)];
          u_chunk[IND(x, y, z - 3)] = u_chunk[IND(x, y, z - 2)];
        }
        if (z == NZr + SHIFT ) {
          u_chunk[IND(x, y, z)] = u_chunk[IND(x, y, z - 1)];
          u_chunk[IND(x, y, z+1)] = u_chunk[IND(x, y, z)];
          u_chunk[IND(x, y, z+2)] = u_chunk[IND(x, y, z + 1)];
          u_chunk[IND(x, y, z+3)] = u_chunk[IND(x, y, z + 2)];
        }
      }
    }
  }
  
  // *************************

 // printf("rank - %d; left %d; right %d; top %d; down %d\n", rankP, rank_left, rank_right, rank_top, rank_down);

  // Создание типа плоскости XY и XZ

  //  TODO:
  //  Подумать над передачей четырёх плокостей!

  //  Каждая плоскость представляет SHIFT подряд идущих плоскостей
  MPI_Datatype planeXY;
  MPI_Type_vector(NZr+RESERVE, NX*SHIFT, NX*(NYr+RESERVE), MPI_TYPE, &planeXY);
  MPI_Type_commit(&planeXY);

  MPI_Datatype planeXZ;
  MPI_Type_contiguous(NX*SHIFT*(NYr+RESERVE), MPI_TYPE, &planeXZ);
  MPI_Type_commit(&planeXZ);
  // *****************************

  if (rankP == ROOT) {
    printf("START!\n");
    t0 = MPI_Wtime();
  }

  // ОСНОВНЫЕ ВЫЧИСЛЕНИЯ
  for (int t = 1; t <= sizeTime; t++) {

    //    Передача влево по Y
    MPI_Sendrecv(&u_chunk[IND(0, NYr, 0)],  1, planeXY, rank_left, 0,
                 &u_chunk[IND(0, 0, 0)], 1, planeXY, rank_right, 0,
                 gridComm, &status[0]);


    //    Передача вправо по Y
    MPI_Sendrecv(&u_chunk[IND(0, SHIFT, 0)], 1, planeXY, rank_right, 1,
                 &u_chunk[IND(0, NYr+SHIFT, 0)],  1, planeXY, rank_left, 1,
                 gridComm, &status[1]);

    //    Передача вниз по Z
    MPI_Sendrecv(&u_chunk[IND(0, 0, 4)], 1, planeXZ, rank_down, 2,
                 &u_chunk[IND(0, 0, NZr+SHIFT)],  1, planeXZ, rank_top, 2,
                 gridComm, &status[2]);

    //    Передача вверх по Z
    MPI_Sendrecv(&u_chunk[IND(0, 0, NZr)], 1, planeXZ, rank_top, 3,
                 &u_chunk[IND(0, 0, 0)], 1, planeXZ, rank_down, 3,
                 gridComm, &status[3]);

    // k1 = A*U
    multMV(k1, A, u_chunk, 0,0,0, coeffs);

    // k2 = B*k1
    multMV(k2, B, k1, 0,0,0, coeffs);

    // k3 = B*k2
    multMV(k3, B, k2, 0,0,0, coeffs);

    // k4 = C*k3
    multMV(k4, C, k3, 0,0,0, coeffs);

    // UNext = U + (k1 + k2*2 + k3*2 + k4)*h;
    sumV(un_chunk, u_chunk, k1, k2, k3, k4, dimChunk, h);

    tmp = u_chunk;
    u_chunk = un_chunk;
    un_chunk = tmp;
  }

  //  *******************
  if (rankP == ROOT) {
    printf("FINISH!\n");
    t1 = MPI_Wtime();
  }

  //  GATHER
  gather_by_block(u, u_chunk, NX, NY, NYr, NZr, RESERVE, gridComm);

  if (rankP == ROOT) {
    double diffTime = t1 - t0;
    printf("Time -\t%.3lf\n", diffTime);
//    writeFunction1D(pathResult, u, NX, NY, 10, 10);
//    writeFunction3D(pathResult3D, u, NX, NY, NZ, SHIFT);
    free(u);
  }
  MPI_Type_free(&planeXY);
  MPI_Type_free(&planeXZ);

  freeSpMat(&A);
  freeSpMat(&B);
  freeSpMat(&C);

  free(un_chunk);
  free(k1);
  free(k2);
  free(k3);
  free(k4);

  MPI_Finalize();
  return 0;
}
