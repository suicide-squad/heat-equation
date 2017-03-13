#include <iostream>
#include <omp.h>
#include "Task.h"

const int ENABLE_PARALLEL = 0;

using std::string;

int main(int argc, char** argv) {

    // Timing variables
    double time_S, time_E;
    string functionFile = "./../../../../../initial/function.txt";
    string settingFile = "./../../../../../initial/setting.ini";

    Task task;
    task.initTaskUsingFile(settingFile, functionFile);

    time_S = omp_get_wtime();
    task.printVectFile("kek", 0);

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


    // TODO тут скорей будет новый алгоритм
//    double expr = (sigma * 1/timesize * (tFinal - tStart)) / (step * step);
//
//    for (double j = 0; j < timesize; j += 1) {
//
//        omp_set_num_threads(4);
//        {
//            #pragma omp parallel for if (ENABLE_PARALLEL)
//            for (int i = 1; i <= nX; i++) {
//                vect[currTime][i] = expr * (vect[prevTime][i + 1] - 2 * vect[prevTime][i] + vect[prevTime][i - 1])
//                                    + vect[prevTime][i];
//            }
//        }
//        // boundary conditions
//        vect[currTime][0] = vect[currTime][1];
//        vect[currTime][nX+1] = vect[currTime][nX];
//
//        prevTime = (prevTime + 1) % 2;
//        currTime = (currTime + 1) % 2;
//
//    }


    time_E = omp_get_wtime();
    printf("Run time:\t %.15lf\n", time_E-time_S);

    // TODO добавить вывод в отдельную функцию
//    string outfilename = "OUTPUT_" + consoleInput + ".txt";
//    FILE *outfile = fopen(outfilename.c_str(), "w");
//
//    for (int i = 1; i <= nX; i++) {
//        fprintf(outfile, "%2.15le\n", vect[prevTime][i]);
//    }
}
