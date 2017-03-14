#include <iostream>
#include <omp.h>
#include "Task.h"
#include "SparseMatrix.h"

using std::string;

int main(int argc, char** argv) {

    // Timing variables
    double time_S, time_E;
    string functionFile = "./../../../../../initial/function.txt";
    string settingFile = "./../../../../../initial/setting.ini";

    Task task;
    task.initTaskUsingFile(settingFile, functionFile);

    time_S = omp_get_wtime();

    TaskExpressions taskexpr;
    task.setTaskExpr(taskexpr);

    task.createMatrix(taskexpr);



    // TODO вынести этот блок в отдельную функцию
//    //printf("%.6lf\n", timedt);
//    string consoleInput = "";
//    if (argv[1] != 0) {
//        timesize = pow(2, atof(argv[1]));
//        consoleInput = argv[1];
//    }

    // TODO а этот ещё в отдельную
//    printf("1/timesize:\t %.10lf\n", 1/timesize);
//    printf("timesize:\t %.0f\n", timesize);
//    printf("dt:\t\t %.2e\n", 1/timesize * (tFinal - tStart));


    // Calculating
    time_S = omp_get_wtime();

    for (double j = 0; j < task.tFinish; j += task.dt) {
        multiplicateVector(task.matrix, task.vect[task.prevTime], task.vect[task.currTime], task.fullVectSize);
        task.prevTime = (task.prevTime + 1) % 2;
        task.currTime = (task.currTime + 1) % 2;
    }
    time_E = omp_get_wtime();
    printf("Run time %.15lf\n", time_E-time_S);


    // Output
    FILE *outfile = fopen("OUTPUT.txt", "w");

    double outData;
    for (int i = 1; i < task.nX; ++i) {
        outData = task.getData(i, 0, 0);
        fprintf(outfile, "%2.15le\n", outData);

    }

//    for (int i = 1; i <= task.fullVectSize; i++) {
//        fprintf(outfile, "%2.15le\n", task.vect[task.prevTime][i]);
//    }
}
