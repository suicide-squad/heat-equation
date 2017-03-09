//
// Created by kirill on 09.03.17.
//

#include <stdlib.h>
#include <stdio.h>

#include "parser.h"

int readSetting(char path[], Setting *setting) {
  FILE *fp;
  if ((fp = fopen(path, "r")) == NULL)
    return NO_FILE;

  if ( !fscanf(fp, "XSTART=%lf\n", &setting->XSTART) )
    return READING_ERROR;
  if ( !fscanf(fp, "XEND=%lf\n", &setting->XEND) )
    return READING_ERROR;
  if ( !fscanf(fp, "YSTART=%lf\n", &setting->YSTART) )
    return READING_ERROR;
  if ( !fscanf(fp, "YEND=%lf\n", &setting->YEND) )
    return READING_ERROR;
  if ( !fscanf(fp, "ZSTART=%lf\n", &setting->ZSTART) )
    return READING_ERROR;
  if ( !fscanf(fp, "ZEND=%lf\n", &setting->ZEND) )
    return READING_ERROR;
  if ( !fscanf(fp, "SIGMA=%lf\n", &setting->SIGMA) )
    return READING_ERROR;
  if ( !fscanf(fp, "NX=%lu\n", &setting->NX) )
    return READING_ERROR;
  if ( !fscanf(fp, "NY=%lu\n", &setting->NY) )
    return READING_ERROR;
  if ( !fscanf(fp, "NZ=%lu\n", &setting->NZ) )
    return READING_ERROR;
  if ( !fscanf(fp, "TSTART=%lf\n", &setting->TSTART) )
    return READING_ERROR;
  if ( !fscanf(fp, "TFINISH=%lf\n", &setting->TFINISH) )
    return READING_ERROR;
  if ( !fscanf(fp, "dt=%lf\n", &setting->dt) )
    return READING_ERROR;
  if ( !fscanf(fp, "BC=%d\n", &setting->BC) )
    return READING_ERROR;
  return OK;
}