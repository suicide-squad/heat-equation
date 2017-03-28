#include <iostream>
#include <omp.h>
#include "Task.h"
#include "SparseMatrix.h"

using std::string;

double getVectorValue(double *vect, int x, int y, int z, Task task);
int main(int argc, char** argv) {

    // Timing variables
    double time_S, time_E;
    int prevTime, currTime;
    string functionFile = "../../initial/function.txt";
    string settingFile = "../../initial/setting.ini";

    double** vect = new double*[2];
    Task task;
    initTaskUsingFile(task, vect, settingFile);

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

    prevTime = 0;
    currTime = 1;

    TaskExpressions taskExpr;
    taskExpr.x1 = (task.sigma * task.dt) / (task.timeStepX * task.timeStepX);
    taskExpr.y1 = (task.sigma * task.dt) / (task.timeStepY * task.timeStepY);
    taskExpr.z1 = (task.sigma * task.dt) / (task.timeStepZ * task.timeStepZ);
    taskExpr.x2Comp = (1 - 2 * taskExpr.x1 - 2 * taskExpr.y1 - 2 * taskExpr.z1);

    SparseMatrix spMat;
    spMatrixInit(spMat, 7*task.nX * task.nY * task.nZ
                        + 2 * task.nY * task.nZ, task.fullVectSize);
    fillMatrix3d6Expr(spMat, taskExpr, task.nX, task.nY, task.nZ);



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
        fprintf(outfile, "%2.15le\n", getVectorValue(vect[0],i,0,0,task));

    }
}

double getVectorValue(double *vect, int x, int y, int z, Task task) {
    return vect[x + (task.nX + 2) * y + (task.nX+2)*task.nY*z];
}

