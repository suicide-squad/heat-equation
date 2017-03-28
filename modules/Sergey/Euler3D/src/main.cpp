#include <iostream>
#include <omp.h>
#include "Task.h"
#include "Task2.h"
#include "SparseMatrix.h"

using std::string;

double getVectorValue(double *vect, int x, int y, int z, Task2 task2);
int main(int argc, char** argv) {

    // Timing variables
    double time_S, time_E;
    int prevTime, currTime;
    string functionFile = "../../initial/function.txt";
    string settingFile = "../../initial/setting.ini";

    double** vect = new double*[2];
    Task2 task;
    initTaskUsingFile(task, vect, settingFile, functionFile);

    FILE *inFunctionfile = fopen(functionFile.c_str(), "r");

    task.fullVectSize = (task.nX + 2) * (task.nY) * (task.nZ);
    vect[0] = new double[task.fullVectSize];
    vect[1] = new double[task.fullVectSize];

    /// Read file
    for (int k = 0; k < task.nZ; k++) {
        for (int j = 0; j < task.nY; ++j) {
            for (int i = 1; i < task.nX + 1; ++i) {
                fscanf(inFunctionfile, "%lf\n", &vect[0][i + (task.nX + 2) * j + (task.nX+2)*task.nY*k]);
            }
        }
    }

    fclose(inFunctionfile);

    preparationData(task);

//    Task task;
//    task.initTaskUsingFile(settingFile, functionFile);
//    task.preparationData();

    prevTime = 0;
    currTime = 1;

    vect[0][0] = vect[0][1];
    vect[0][task.nX+1] = vect[0][task.nX];

    TaskExpressions taskexpr;
    setTaskExpr(taskexpr, task);

    SparseMatrix spMat;
    spMatrixInit(spMat, 7*task.nX * task.nY * task.nZ
                        + 2 * task.nY * task.nZ, task.fullVectSize);
    fillMatrix3d6Expr(spMat, taskexpr, task.nX, task.nY, task.nZ);
//    task.setTaskExpr(taskexpr);
//    task.createMatrix(taskexpr);



    // Calculating
    time_S = omp_get_wtime();

    for (double j = 0; j < task.tFinish; j += task.dt) {
//        multiplicateVector(task.matrix, task.vect[task.prevTime], task.vect[task.currTime], task.fullVectSize);
        multiplicateVector(spMat, vect[prevTime], vect[currTime], task.fullVectSize);
        prevTime = (prevTime + 1) % 2;
        currTime = (currTime + 1) % 2;
    }
    time_E = omp_get_wtime();
    printf("Run time %.15lf\n", time_E-time_S);


    // Output
    FILE *outfile = fopen("../../result/Sergey/Sergey_Euler_OUTPUT.txt", "w");

    double outData;
    for (int i = 1; i <= task.nX; ++i) {
        fprintf(outfile, "%2.15le\n", getVectorValue(vect[0],0,0,0,task));

    }
}

double getVectorValue(double *vect, int x, int y, int z, Task2 task2) {
    return vect[x + (task2.nX + 2) * y + (task2.nX+2)*task2.nY*z];
}

