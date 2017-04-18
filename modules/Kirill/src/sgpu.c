//
// Created by kirill on 12.04.17.
//

#include <stdlib.h>
#include <math.h>

#include "sgpu.h"

void scatter_by_block(double *u, double *u_chunk, int NX, int NY, int NYr, int NZr, MPI_Comm gridComm) {
  int ndim = 0;
  MPI_Cartdim_get(gridComm, &ndim);
  int dims[ndim], periods[ndim], gridCoords[ndim];

  MPI_Cart_get(gridComm, ndim, dims, periods, gridCoords);
  MPI_Comm colComm, rowComm;

  int remain_dims[ndim];
  remain_dims[0] = 1; remain_dims[1] = 0;
  MPI_Cart_sub(gridComm, remain_dims, &colComm);

  double *u_z = NULL;

  if (gridCoords[1] == 0) {
    u_z = (double *)malloc(sizeof(double)*NX*NY*(NZr+2));
    MPI_Scatter(u + NY*NX, NX*NY*(NZr), MPI_DOUBLE, u_z + NX*NY, NX*NY*(NZr), MPI_DOUBLE, 0, colComm);
  }

  remain_dims[0] = 0; remain_dims[1] = 1;
  MPI_Cart_sub(gridComm, remain_dims, &rowComm);

  for (int i = 1; i < NZr + 1; i++)
    MPI_Scatter(u_z + NX + i*NX*NY, NX*(NYr), MPI_DOUBLE, u_chunk + NX + i*NX*(NYr+2), NX*NYr, MPI_DOUBLE, 0, rowComm);

  if (gridCoords[1] == 0) free(u_z);
}


void gather_by_block(double *u, double *u_chunk, int NX, int NY, int NYr, int NZr, MPI_Comm gridComm) {
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

  double *u_z = NULL;
  if (gridCoords[1] == 0) u_z = (double *) malloc(sizeof(double) * NX * NY *(NZr));

  for (int z = 0; z < NZr; z++)
    MPI_Gather(u_chunk + NX + (z+1)*NX*(NYr+2), NX*(NYr), MPI_DOUBLE, u_z + NX + z*NX*NY, NX*(NYr), MPI_DOUBLE, 0, rowComm);

  if (gridCoords[1] == 0) {
    MPI_Gather(u_z, NX*NY*(NZr), MPI_DOUBLE, u + NY*NX, NX*NY*(NZr), MPI_DOUBLE, 0, colComm);
    free(u_z);
  }
}
void get_blocks(int* blockY, int* blockZ, int sizeProc) {
  int sqrtNum = (int)sqrt(sizeProc);
  for(int i = sqrtNum; i <= sizeProc; i++) {
    if (sizeProc % i == 0) {
      *blockY = i;
      break;
    }
  }
  *blockZ = sizeProc/(*blockY);
}
