//
// Created by kirill on 12.04.17.
//

#include <stdlib.h>

#include "sgpu.h"

void pack(double *ex, double *u, int NX, int NY, int NZ, op which) {
  switch ( which ) {
    case Y_LEFT_SEND: {
      for (int z = 0; z < NZ; z++)
        for (int x = 0; x < NX; x++)
          ex[x + z * NX] = u[x + 1 * NX + z * NX * NZ];
      break;
    }

    case Y_RIGHT_SEND: {
      for (int z = 0; z < NZ; z++)
        for (int x = 0; x < NX; x++)
          ex[x + z * NX] = u[x + (NY - 2) * NX + z * NX * NZ];
      break;
    }

    case Z_DOWN_SEND: {
      for (int y = 0; y < NY; y++)
        for (int x = 0; x < NX; x++)
          ex[x + y * NX] = u[x + y * NX + 1 * NX * NY];
      break;
    }
    case Z_TOP_SEND: {
      for (int y = 0; y < NY; y++)
        for (int x = 0; x < NX; x++)
          ex[x + y * NX] = u[x + y * NX + (NZ - 2) * NX * NY];
      break;
    }
  }
  return;
}

void unpack (double *ex, double *u, int NX, int NY, int NZ, op which) {
  switch ( which ) {
    case Y_LEFT_RECV: {
      for (int z = 0; z < NZ; z++)
        for (int x = 0; x < NX; x++)
          u[x + 0 * NX + z * NX * NZ] = ex[x + z * NX];
      break;
    }
    case Y_RIGHT_RECV: {
      for (int z = 0; z < NZ; z++)
        for (int x = 0; x < NX; x++)
          u[x + (NY - 1) * NX + z * NX * NZ] = ex[x + z * NX];
      break;
    }
    case Z_DOWN_RECV: {
      for (int y = 0; y < NY; y++)
        for (int x = 0; x < NX; x++)
          u[x + y * NX + 0 * NX * NY] = ex[x + y * NX];
      break;
    }
    case Z_TOP_RECV: {
      for (int y = 0; y < NY; y++)
        for (int x = 0; x < NX; x++)
          u[x + y * NX + (NZ - 1) * NX * NY] = ex[x + y * NX];
      break;
    }
  }
  return;
}


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
    MPI_Scatter(u, NX*NY*(NZr), MPI_DOUBLE, u_z + NX*NY, NX*NY*(NZr), MPI_DOUBLE, 0, colComm);
  }

  remain_dims[0] = 0; remain_dims[1] = 1;
  MPI_Cart_sub(gridComm, remain_dims, &rowComm);

  for (int i = 1; i < NZr + 1; i++)
    MPI_Scatter(u_z + i*NX*NY, NX*(NYr), MPI_DOUBLE, u_chunk + NX + i*NX*(NYr), NX*NYr, MPI_DOUBLE, 0, rowComm);

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
  if (gridCoords[1] == 0) u_z = (double *) malloc(sizeof(double) * NX * NY *NZr);

  for (int i = 0; i < NZr; i++)
    MPI_Gather(u_chunk + NX + (i+1)*(NYr)*NX, NX*(NYr), MPI_DOUBLE, u_z + i*NX*NY, NX*(NYr), MPI_DOUBLE, 0, rowComm);

  if (gridCoords[1] == 0) {
    MPI_Gather(u_z, NX * NY * (NZr), MPI_DOUBLE,
               u, NX * NY * (NZr), MPI_DOUBLE, 0,
               colComm);
    free(u_z);
  }
  int rank;
  MPI_Comm_rank(gridComm, &rank);
}