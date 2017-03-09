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
  if ( !fscanf(fp, "NX=%lu\n", &setting->NX) ) return READING_ERROR;
  if ( !fscanf(fp, "NY=%lu\n", &setting->NY) ) return READING_ERROR;
  if ( !fscanf(fp, "NZ=%lu\n", &setting->NZ) ) return READING_ERROR;
  if ( !fscanf(fp, "TSTART=%lf\n", &setting->TSTART) ) return READING_ERROR;
  if ( !fscanf(fp, "TFINISH=%lf\n", &setting->TFINISH) ) return READING_ERROR;
  if ( !fscanf(fp, "dt=%lf\n", &setting->dt) ) return READING_ERROR;
  if ( !fscanf(fp, "BC=%d\n", &setting->BC) ) return READING_ERROR;
  return OK;
}

int readFunction(const char *path, double **function, size_t dim, size_t NX) {
  FILE *fp;
  if ( !(fp = fopen(path, "r")) )
    return NO_FILE;

  for (int i = 0; i < dim; i++)
    if ( !fscanf(fp, "%lf\n", (*function) + i)
        && i%(NX + 2)!=0 && i%(NX + 2) != NX + 1 )
      return READING_ERROR;

  return OK;
}

int writeFunctionX(const char *path, double *function, size_t dim, size_t NX) {
  FILE *fp;
  fp = fopen(path, "w");

  for (int i = 1; i < NX + 1; i++)
    fprintf(fp, "%.15le\n", function[i]);

  fclose(fp);
  return OK;
}
