//
// Created by kirill on 09.03.17.
//

#ifndef HEAT_EQUATION_PARSER_H
#define HEAT_EQUATION_PARSER_H

#include "utils/ts.h"

#define OK            0
#define NO_FILE       100
#define READING_ERROR 101

typedef struct {
  TYPE XSTART;
  TYPE XEND;
  TYPE YSTART;
  TYPE YEND;
  TYPE ZSTART;
  TYPE ZEND;
  TYPE SIGMA;
  int NX;
  int NY;
  int NZ;
  TYPE TSTART;
  TYPE TFINISH;
  TYPE dt;
  int BC;
} Setting;

int readSetting(const char *path, Setting *setting);
int readFunction(const char *path, TYPE *u, int nx, int ny, int nz, int shift);
void writeFunction1D(const char *path, TYPE *u, int nx, int ny, int y, int z);
void writeFunction3D(const char *path, TYPE *u, int nx, int ny, int nz, int shift);

#endif //HEAT_EQUATION_PARSER_H
