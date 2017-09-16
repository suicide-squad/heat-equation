//
// Created by kirill on 17.09.17.
//

#ifndef HEAT_EQUATION_TS_H_H
#define HEAT_EQUATION_TS_H_H

#include <mpi.h>

#define ENABLE_PARALLEL 0
#define ROOT 0

//TODO: Как грамотно считать тип?
#define DOUBLE_TYPE

#if defined(COMPLEX_TYPE)

#include <complex.h>
#define DEF_TYPE COMPLEX
typedef complex TYPE;

#elif defined(FLOAT_TYPE)

#define MPI_TYPE MPI_FLOAT
typedef float TYPE;

#elif defined(DOUBLE_TYPE)

#define MPI_TYPE MPI_DOUBLE
typedef double TYPE;

#endif

#endif //HEAT_EQUATION_TS_H_H
