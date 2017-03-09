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
    int initTaskUsingFile(string filename);
    int initMemory(int size = nX);
    void printTaskData();
    void preparationData();

    int     sizetime;
    int     currTime, prevTime;

private:
    double** vect;
    double  xStart, xEnd;   // range
    double  sigma;          //
    int     bc;             // not used

    int     nX;             // count of initial elements
    double  tStart, tFinish;
    double  dt;

    double  stepSize;       // distance between two dot
};


#endif //EULER_TASK_H