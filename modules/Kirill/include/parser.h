//
// Created by kirill on 09.03.17.
//

#ifndef HEAT_EQUATION_PARSER_H
#define HEAT_EQUATION_PARSER_H

#define OK 0
#define NO_FILE -100
#define READING_ERROR -101


typedef struct {
  double XSTART;
  double XEND;
  double YSTART;
  double YEND;
  double ZSTART;
  double ZEND;
  double SIGMA;
  size_t NX;
  size_t NY;
  size_t NZ;
  double TSTART;
  double TFINISH;
  double dt;
  int BC;
} Setting;

int readSetting(char *path, Setting *setting);
void readFunction(char path[], Setting* setting);

#endif //HEAT_EQUATION_PARSER_H
