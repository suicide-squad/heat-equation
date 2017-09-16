//
// Created by kirill on 12.04.17.
//

#include <stdlib.h>
#include <math.h>

#include "sgpu.h"

void scatter_by_block(TYPE* u, TYPE* u_chunk, int nx, int ny, int nyr, int nzr, MPI_Comm gridComm, int reserve) {
  int ndim = 0;
  MPI_Cartdim_get(gridComm, &ndim);
  int dims[ndim], periods[ndim], gridCoords[ndim];

  int indRsrv = reserve/2;

  MPI_Cart_get(gridComm, ndim, dims, periods, gridCoords);
  MPI_Comm colComm, rowComm;

  int remain_dims[ndim];
  remain_dims[0] = 1; remain_dims[1] = 0;
  MPI_Cart_sub(gridComm, remain_dims, &colComm);

  TYPE *u_z = NULL;

  if (gridCoords[1] == 0) {
    u_z = (TYPE *)malloc(sizeof(TYPE)*nx*ny*(nzr+reserve));
    MPI_Scatter(u + indRsrv*ny*nx, nx*ny*(nzr), MPI_TYPE, u_z + indRsrv*nx*ny, nx*ny*(nzr), MPI_TYPE, 0, colComm);
  }

  remain_dims[0] = 0; remain_dims[1] = 1;
  MPI_Cart_sub(gridComm, remain_dims, &rowComm);

  for (int i = indRsrv; i < nzr + indRsrv; i++)
    MPI_Scatter(u_z + indRsrv*nx + i*nx*ny, nx*(nyr), MPI_TYPE, u_chunk + indRsrv*nx + i*nx*(nyr+reserve), nx*nyr, MPI_TYPE, 0, rowComm);

  if (gridCoords[1] == 0) free(u_z);
}


void gather_by_block(TYPE* u, TYPE* u_chunk, int nx, int ny, int nyr, int nzr, int reserve, MPI_Comm gridComm) {
  int indRsrv = reserve/2;

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

  TYPE *u_z = NULL;
  if (gridCoords[1] == 0) u_z = (TYPE *)malloc(sizeof(TYPE)*nx*ny*(nzr));

  for (int z = 0; z < nzr; z++)
    MPI_Gather(u_chunk + indRsrv*nx + (z+indRsrv)*nx*(nyr+reserve),
               nx*(nyr), MPI_TYPE, u_z + indRsrv*nx + z*nx*ny, nx*(nyr), MPI_TYPE, 0, rowComm);

  if (gridCoords[1] == 0) {
    MPI_Gather(u_z, nx*ny*(nzr), MPI_TYPE, u + indRsrv*ny*nx, nx*ny*(nzr), MPI_TYPE, 0, colComm);
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
