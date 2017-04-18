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

    // File variables
//    string functionFile = "../../initial/function.txt";
    string functionFile = "../../initial/function2.txt";
//    string functionFile = "../../initial_test/function.txt";
//    string settingFile = "../../initial/setting.ini";
    string settingFile = "../../initial/setting2.ini";
//    string settingFile = "../../initial_test/setting.ini";

    // Read task settings
    Task task;
    initTaskUsingFile(task, settingFile);
    setTimestep(task);

    // Init memory & read function file
    double** vect;
    initMemoryReadData(vect, functionFile, task);

    // vector time-index for loop
    prevTime = 0;
    currTime = 1;

    // value for the matrix
    MatrixValue matrixValue;
    matrixValue.x1 = (task.sigma * task.dt) / (task.timeStepX * task.timeStepX);
    matrixValue.y1 = (task.sigma * task.dt) / (task.timeStepY * task.timeStepY);
    matrixValue.z1 = (task.sigma * task.dt) / (task.timeStepZ * task.timeStepZ);
    matrixValue.x2Comp = (1 - 2 * matrixValue.x1 - 2 * matrixValue.y1 - 2 * matrixValue.z1);

    // init and fill sparseMatrix
    SparseMatrix spMat;
    int sparseMatrixSize = 9 * task.nX * task.nY * task.nZ;

    spMatrixInit(spMat, sparseMatrixSize, task.fullVectSize);
    fillMatrix3d6Expr(spMat, matrixValue, task.nX, task.nY, task.nZ);



    // Calculating
    time_S = omp_get_wtime();

    for (double j = 0; j < task.tFinish; j += task.dt) {
        multiplicateVector(spMat, vect[prevTime], vect[currTime], task.fullVectSize);
        prevTime = (prevTime + 1) % 2;
        currTime = (currTime + 1) % 2;
    }
    time_E = omp_get_wtime();
    printf("Run time %.15lf\n", time_E-time_S);


    // Output
    FILE *outfile = fopen("../../result/Sergey/Sergey_Euler2_full.txt", "w");

    for (int i = 0; i < task.fullVectSize; ++i) {
        if ( i % (task.nX + 2) != 0 && i % (task.nX + 2) != task.nX + 1)
            fprintf(outfile, "%2.15le\n", vect[0][i]);
    }
//    double outData;
//    for (int i = 1; i <= task.nX; ++i) {
//        fprintf(outfile, "%2.15le\n", getVectorValue(vect[0],i,0,0,task));
//
//    }
}

double getVectorValue(double *vect, int x, int y, int z, Task task) {
    return vect[x + (task.nX + 2) * y + (task.nX+2)*task.nY*z];
}

