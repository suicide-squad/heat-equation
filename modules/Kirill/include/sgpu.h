//
// Created by kirill on 12.04.17.
//

#ifndef HEAT_EQUATION_SGPU_H
#define HEAT_EQUATION_SGPU_H

#include <mpi.h>

typedef enum {
  Y_RIGHT_SEND,
  Y_RIGHT_RECV,
  Y_LEFT_SEND,
  Y_LEFT_RECV,
  Z_TOP_SEND,
  Z_TOP_RECV,
  Z_DOWN_SEND,
  Z_DOWN_RECV
} op;

void get_blocks(int *blockY, int *blockZ, int sizeProc);

void scatter_by_block(double *u, double *u_chunk, int NX, int NY, int NYr, int NZr, MPI_Comm gridComm);
void gather_by_block(double *u, double *u_chunk, int NX, int NY, int NYr, int NZr, MPI_Comm gridComm);


#endif //HEAT_EQUATION_SGPU_H
