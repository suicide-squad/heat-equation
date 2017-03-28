//
// Created by lenferd on 27.03.17.
//

#ifndef HEAT_EQUATION_TASK_H
#define HEAT_EQUATION_TASK_H

#include <iostream>
using std::string;

struct Task {
    double  xStart, xEnd;   // range

    double  yStart, yEnd;
    double  zStart, zEnd;

    double  sigma;          //
    int     bc;             // not used

    double  timeStepX;       // time time between calculating
    double  timeStepY;
    double  timeStepZ;

    int     nX;             // count of initial elements
    int     nY;
    int     nZ;
    int     fullVectSize;

    double  tStart, tFinish;
    double  dt;
};

int initTaskUsingFile(Task &task, double **vect, string settingFile);
int preparationData(Task &task);
#endif //HEAT_EQUATION_TASK_H
