//
// Created by kirill on 09.03.17.
//

#ifndef HEAT_EQUATION_PARSER_H
#define HEAT_EQUATION_PARSER_H

#define OK 0
#define NO_FILE 100
#define READING_ERROR 101

typedef struct {
  double XSTART;
  double XEND;
  double YSTART;
  double YEND;
  double ZSTART;
  double ZEND;
  double SIGMA;
  int NX;
  int NY;
  int NZ;
  double TSTART;
  double TFINISH;
  double dt;
  int BC;
} Setting;

int readSetting(const char *path, Setting *setting);
int readFunction(const char* path, double* u, int nx, int ny, int nz, int shift);
void writeFunction1D(const char *path, double *u, int nx, int ny, int y, int z);
void writeFunction3D(const char* path, double* u, int nx, int ny, int nz, int shift);

#endif //HEAT_EQUATION_PARSER_H
