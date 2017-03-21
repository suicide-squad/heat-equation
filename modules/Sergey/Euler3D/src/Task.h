//
// Created by lenferd on 09.03.17.
//

#ifndef EULER_TASK_H
#define EULER_TASK_H

#include "SparseMatrix.h"
#include "StructDeclamer.h"
#include <string>
#include <cmath>

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

    double  tStart, tFinish;
    double  dt;
    int     currTime, prevTime;

    int fullVectSize;
    int     nX;             // count of initial elements
    int     nY;
    int     nZ;

    // Sparse Matrix
    SparseMatrix matrix;

    double** vect;

private:

    double  xStart, xEnd;   // range

    double  yStart, yEnd;
    double  zStart, zEnd;

    double  sigma;          //
    int     bc;             // not used
    double  timeStepX;       // time time between calculating

    double  timeStepY;

    double  timeStepZ;

};

#endif //EULER_TASK_H
