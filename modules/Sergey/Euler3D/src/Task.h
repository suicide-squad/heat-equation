//
// Created by lenferd on 09.03.17.
//

#ifndef EULER_TASK_H
#define EULER_TASK_H

#include <string>
#include <cmath>

using std::string;

class Task {
public:
    Task();
    int initTaskUsingFile(string settingFile, string functionFile);
    int initMemory();
    void preparationData();
    void setData(double data, int x, int y, int z, int timeVect = 0);

    double getData(int x, int y, int z, int timeVect = 0);

    void printTaskData();
    void printVect(int timeVect = 0);

    // debug function
    void printVectFile(string filename, int timeVect = 0);

    int     sizetime;
    int     currTime, prevTime;

private:
    double** vect;
    double  xStart, xEnd;   // range
    double  yStart, yEnd;
    double  zStart, zEnd;
    double  sigma;          //
    int     bc;             // not used

    int     nX;             // count of initial elements
    int     nY;
    int     nZ;

    double  tStart, tFinish;
    double  dt;

    double  stepSize;       // distance between two dot

    int fullVectSize;
};


#endif //EULER_TASK_H
