//
// Created by kirill on 09.03.17.
//

#include <stdlib.h>
#include <stdio.h>

#include "parser.h"

int readSetting(const char *path, Setting *setting) {
  FILE *fp;
  if ( !(fp = fopen(path, "r")) )
    return NO_FILE;

  if ( !fscanf(fp, "XSTART=%lf\n", &setting->XSTART) ) return READING_ERROR;
  if ( !fscanf(fp, "XEND=%lf\n", &setting->XEND) ) return READING_ERROR;
  if ( !fscanf(fp, "YSTART=%lf\n", &setting->YSTART) ) return READING_ERROR;
  if ( !fscanf(fp, "YEND=%lf\n", &setting->YEND) ) return READING_ERROR;
  if ( !fscanf(fp, "ZSTART=%lf\n", &setting->ZSTART) ) return READING_ERROR;
  if ( !fscanf(fp, "ZEND=%lf\n", &setting->ZEND) ) return READING_ERROR;
  if ( !fscanf(fp, "SIGMA=%lf\n", &setting->SIGMA) ) return READING_ERROR;
  if ( !fscanf(fp, "NX=%d\n", &setting->NX) ) return READING_ERROR;
  if ( !fscanf(fp, "NY=%d\n", &setting->NY) ) return READING_ERROR;
  if ( !fscanf(fp, "NZ=%d\n", &setting->NZ) ) return READING_ERROR;
  if ( !fscanf(fp, "TSTART=%lf\n", &setting->TSTART) ) return READING_ERROR;
  if ( !fscanf(fp, "TFINISH=%lf\n", &setting->TFINISH) ) return READING_ERROR;
  if ( !fscanf(fp, "dt=%lf\n", &setting->dt) ) return READING_ERROR;
  if ( !fscanf(fp, "BC=%d\n", &setting->BC) ) return READING_ERROR;
  return OK;
}

int readFunction(const char* path, double* u, int nx, int ny, int nz, int shift) {
  FILE *fp;
  if ( !(fp = fopen(path, "r")) )
    return NO_FILE;

  for (int z = shift; z < nz - shift; z++)
    for (int y = shift; y < ny - shift; y++)
      for (int x = shift; x < nx - shift; x++)
        if (!fscanf(fp, "%lf\n", &u[x + y*nx + z*nx*ny]))
          return READING_ERROR;

  return OK;
}

void writeFunction1D(const char *path, double *u, int nx, int ny, int y, int z) {
  FILE *fp;
  fp = fopen(path, "w");

  for (int i = 1; i < nx - 1; i++)
    fprintf(fp, "%.15le\n", u[i+nx*y+nx*ny*z]);

  fclose(fp);
}

void writeFunction3D(const char* path, double* u, int nx, int ny, int nz, int shift) {
  FILE *fp;
  fp = fopen(path, "w");

  for (int z = shift; z < nz-shift; z++)
    for (int y = shift; y < ny-shift; y++)
      for (int x = shift; x < nx-shift; x++)
        fprintf(fp, "%.15le\n", u[x + y*nx + z*nx*ny]);
  fclose(fp);
}
