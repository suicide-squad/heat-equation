//
// Created by lenferd on 28.03.17.
//

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
    string functionFile = "../../initial/function.txt";
    string settingFile = "../../initial/setting.ini";

    // Read task settings
    Task task;
    initTaskUsingFile(task, settingFile);
    setTimestep(task);

    // Init memory & read function file
    double** vect;
    initMemoryReadData(vect, functionFile, task);


    double* vectK1 = new double[task.fullVectSize];
    double* vectK2 = new double[task.fullVectSize];
    double* vectK3 = new double[task.fullVectSize];
    double* vectK4 = new double[task.fullVectSize];


    // vector time-index for loop
    prevTime = 0;
    currTime = 1;

    /***
     * K1 Vector
     */
    // value for the matrix
    MatrixValue matrixValueK1;
    matrixValueK1.x1 = (task.sigma) / (task.timeStepX * task.timeStepX);
    matrixValueK1.y1 = (task.sigma) / (task.timeStepY * task.timeStepY);
    matrixValueK1.z1 = (task.sigma) / (task.timeStepZ * task.timeStepZ);
    matrixValueK1.x2Comp = (- 2 * matrixValueK1.x1 - 2 * matrixValueK1.y1 - 2 * matrixValueK1.z1);

    // init and fill sparseMatrix
    SparseMatrix spMatK1;
    int sparseMatrixSize = 9 * task.nX * task.nY * task.nZ;

    spMatrixInit(spMatK1, sparseMatrixSize, task.fullVectSize);
    fillMatrix3d6Expr(spMatK1, matrixValueK1, task.nX, task.nY, task.nZ);

    /***
    * K2-3 Vector
    */
    // value for the matrix
    MatrixValue matrixValueK2;
    matrixValueK2.x1 = (task.sigma * task.dt * 0.5) / (task.timeStepX * task.timeStepX);
    matrixValueK2.y1 = (task.sigma * task.dt * 0.5) / (task.timeStepY * task.timeStepY);
    matrixValueK2.z1 = (task.sigma * task.dt * 0.5) / (task.timeStepZ * task.timeStepZ);
    matrixValueK2.x2Comp = (1 - 2 * matrixValueK2.x1 - 2 * matrixValueK2.y1 - 2 * matrixValueK2.z1);

    // init and fill sparseMatrix
    SparseMatrix spMatK2;
    sparseMatrixSize = 9 * task.nX * task.nY * task.nZ;

    spMatrixInit(spMatK2, sparseMatrixSize, task.fullVectSize);
    fillMatrix3d6Expr(spMatK2, matrixValueK2, task.nX, task.nY, task.nZ);

    /***
    * K4 Vector
    */
    // value for the matrix
    MatrixValue matrixValueK4;
    matrixValueK4.x1 = (task.sigma * task.dt) / (task.timeStepX * task.timeStepX);
    matrixValueK4.y1 = (task.sigma * task.dt) / (task.timeStepY * task.timeStepY);
    matrixValueK4.z1 = (task.sigma * task.dt) / (task.timeStepZ * task.timeStepZ);
    matrixValueK4.x2Comp = (1 - 2 * matrixValueK4.x1 - 2 * matrixValueK4.y1 - 2 * matrixValueK4.z1);

    // init and fill sparseMatrix
    SparseMatrix spMatK4;
    sparseMatrixSize = 9 * task.nX * task.nY * task.nZ;

    spMatrixInit(spMatK4, sparseMatrixSize, task.fullVectSize);
    fillMatrix3d6Expr(spMatK4, matrixValueK4, task.nX, task.nY, task.nZ);




    // Calculating
    time_S = omp_get_wtime();

    for (double j = 0; j < task.tFinish; j += task.dt) {
        multiplicateVector(spMatK1, vect[prevTime], vectK1, task.fullVectSize);
        multiplicateVector(spMatK2, vectK1, vectK2, task.fullVectSize);
        multiplicateVector(spMatK2, vectK2, vectK3, task.fullVectSize); // Sparse matrix for k3 equal k2
        multiplicateVector(spMatK4, vectK3, vectK4, task.fullVectSize);

        for (int i = 0; i < task.fullVectSize; ++i) {
            vect[currTime][i] =
                    vect[prevTime][i] + task.dt / 6 * (vectK1[i] + 2 * vectK2[i] + 2 * vectK3[i] + vectK4[i]);
        }

        prevTime = (prevTime + 1) % 2;
        currTime = (currTime + 1) % 2;
    }
    time_E = omp_get_wtime();
    printf("Run time %.15lf\n", time_E-time_S);


    // Output
    FILE *outfile = fopen("../../result/Sergey/Sergey_Runge.txt", "w");

    double outData;
    for (int i = 1; i <= task.nX; ++i) {
        fprintf(outfile, "%2.15le\n", getVectorValue(vect[0],i,0,0,task));

    }
}

double getVectorValue(double *vect, int x, int y, int z, Task task) {
    return vect[x + (task.nX + 2) * y + (task.nX+2)*task.nY*z];
}

