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

  if ( !fscanf(fp, "XSTART=%" PRI_TYPE "\n", &setting->XSTART) ) return READING_ERROR;
  if ( !fscanf(fp, "XEND=%" PRI_TYPE "\n", &setting->XEND) ) return READING_ERROR;
  if ( !fscanf(fp, "YSTART=%" PRI_TYPE "\n", &setting->YSTART) ) return READING_ERROR;
  if ( !fscanf(fp, "YEND=%" PRI_TYPE "\n", &setting->YEND) ) return READING_ERROR;
  if ( !fscanf(fp, "ZSTART=%" PRI_TYPE "\n", &setting->ZSTART) ) return READING_ERROR;
  if ( !fscanf(fp, "ZEND=%" PRI_TYPE "\n", &setting->ZEND) ) return READING_ERROR;
  if ( !fscanf(fp, "SIGMA=%" PRI_TYPE "\n", &setting->SIGMA) ) return READING_ERROR;
  if ( !fscanf(fp, "NX=%d\n", &setting->NX) ) return READING_ERROR;
  if ( !fscanf(fp, "NY=%d\n", &setting->NY) ) return READING_ERROR;
  if ( !fscanf(fp, "NZ=%d\n", &setting->NZ) ) return READING_ERROR;
  if ( !fscanf(fp, "TSTART=%" PRI_TYPE "\n", &setting->TSTART) ) return READING_ERROR;
  if ( !fscanf(fp, "TFINISH=%" PRI_TYPE "\n", &setting->TFINISH) ) return READING_ERROR;
  if ( !fscanf(fp, "dt=%" PRI_TYPE "\n", &setting->dt) ) return READING_ERROR;
  if ( !fscanf(fp, "BC=%d\n", &setting->BC) ) return READING_ERROR;
  return OK;
}

int readFunction(const char *path, TYPE *u, int nx, int ny, int nz, int shift) {
  FILE *fp;
  if ( !(fp = fopen(path, "r")) )
    return NO_FILE;

  for (int z = shift; z < nz - shift; z++)
    for (int y = shift; y < ny - shift; y++)
      for (int x = shift; x < nx - shift; x++)
        if (!fscanf(fp, "%" PRI_TYPE "\n", &u[x + y*nx + z*nx*ny]))
          return READING_ERROR;

  return OK;
}

void writeFunction1D(const char *path, TYPE *u, int nx, int ny, int y, int z) {
  FILE *fp;
  fp = fopen(path, "w");

  for (int i = 1; i < nx - 1; i++)
    fprintf(fp, "%" PRIE_TYPE "\n", u[i+nx*y+nx*ny*z]);

  fclose(fp);
}

void writeFunction3D(const char *path, TYPE *u, int nx, int ny, int nz, int shift) {
  FILE *fp;
  fp = fopen(path, "w");

  for (int z = shift; z < nz-shift; z++)
    for (int y = shift; y < ny-shift; y++)
      for (int x = shift; x < nx-shift; x++)
        fprintf(fp, "%" PRIE_TYPE "\n", u[x + y*nx + z*nx*ny]);
  fclose(fp);
}
