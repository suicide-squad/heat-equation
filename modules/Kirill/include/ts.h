//
// Created by kirill on 17.09.17.
//

#ifndef HEAT_EQUATION_TS_H_H
#define HEAT_EQUATION_TS_H_H

#include <mpi.h>
#include <math.h>

#define ENABLE_PARALLEL 0
#define ROOT 0

//TODO: Как грамотно считать тип?
#define FLOAT_TYPE

#if defined(COMPLEX_TYPE)

#include <complex.h>
#define DEF_TYPE COMPLEX
typedef complex TYPE;

#elif defined(FLOAT_TYPE)

#define MPI_TYPE MPI_FLOAT
#define PRI_TYPE "f"
#define PRIE_TYPE ".7e"
#define ABS(x) fabsf(x)
typedef float TYPE;

#elif defined(DOUBLE_TYPE)

#define MPI_TYPE MPI_DOUBLE
#define PRI_TYPE "lf"
#define PRIE_TYPE ".15le"
#define ABS(x) fabs(x)
typedef double TYPE;

#endif

#endif //HEAT_EQUATION_TS_H_H
