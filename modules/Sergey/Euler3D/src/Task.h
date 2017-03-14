//
// Created by lenferd on 09.03.17.
//

#ifndef EULER_TASK_H
#define EULER_TASK_H

#include <string>
#include <cmath>
#include "SparseMatrix.h"

using std::string;

class Task {
public:
    Task();

    // first step
    int initTaskUsingFile(string settingFile, string functionFile);
    int initMemory();
    void preparationData();

    // get and set
    void setData(double data, int x, int y, int z, int timeVect = 0);
    double getData(int x, int y, int z, int timeVect = 0);

    // Transform
    void createMatrix(TaskExpressions &taskexpr);
    void setTaskExpr(TaskExpressions &task);

    // Print function
    void printTaskData();
    void printVect(int timeVect = 0);

    // debug function
    void printVectFile(string filename, int timeVect = 0);

    int     currTime, prevTime;

    int     nX;             // count of initial elements
    int     nY;
    int     nZ;

private:
    double** vect;
    double  xStart, xEnd;   // range
    double  yStart, yEnd;
    double  zStart, zEnd;
    double  sigma;          //

    int     bc;             // not used

    double  tStart, tFinish;
    double  dt;

    double  timeStepX;       // time time between calculating
    double  timeStepY;
    double  timeStepZ;

    int fullVectSize;

    // Sparse Matrix
    SparseMatrix matrix;

};

typedef struct TaskExpressions {
    double x1;
    double x2Comp;

    double y1;
//    double y2;

    double z1;
//    double z2;
};

#endif //EULER_TASK_H
