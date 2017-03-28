//
// Created by lenferd on 27.03.17.
//

#include <cmath>
#include "Task2.h"

using std::string;

int initTaskUsingFile(Task2 &task, double **vect, string settingFile, string functionFile) {
    FILE *inSettingfile = fopen(settingFile.c_str(), "r");

    if (inSettingfile == NULL) {
        printf("File reading error. Try to relocate input file\n");
        exit(0);
    }

//    XSTART=-1.0
//    XEND=1.0
//    YSTART=-1.0
//    YEND=1.0
//    ZSTART=-1.0
//    ZEND=1.0
//    SIGMA=1.0
//    NX=50
//    NY=50
//    NZ=50
//    TSTART=0.000
//    TFINISH=8e-3
//    dt=8e-8
//    BC=2

    // File reading
    fscanf(inSettingfile, "XSTART=%lf\n", &task.xStart);    // start coordinate
    fscanf(inSettingfile, "XEND=%lf\n", &task.xEnd);        // end coordinate

    fscanf(inSettingfile, "YSTART=%lf\n", &task.yStart);    // start coordinate
    fscanf(inSettingfile, "YEND=%lf\n", &task.yEnd);        // end coordinate

    fscanf(inSettingfile, "ZSTART=%lf\n", &task.zStart);    // start coordinate
    fscanf(inSettingfile, "ZEND=%lf\n", &task.zEnd);        // end coordinate

    fscanf(inSettingfile, "SIGMA=%lf\n", &task.sigma);      // coef of heat conduction

    fscanf(inSettingfile, "NX=%d\n", &task.nX);             // count of initial elements
    fscanf(inSettingfile, "NY=%d\n", &task.nY);             //
    fscanf(inSettingfile, "NZ=%d\n", &task.nZ);             //

    fscanf(inSettingfile, "TSTART=%lf\n", &task.tStart);    // start time
    fscanf(inSettingfile, "TFINISH=%lf\n", &task.tFinish);   // finish time
    fscanf(inSettingfile, "dt=%lf\n", &task.dt);            // delta of time difference
    fscanf(inSettingfile, "BC=%d\n", &task.bc);         // Not using right now

    fclose(inSettingfile);
    return 0;
}

int preparationData(Task2 &task){
    task.timeStepX = (fabs(task.xStart) + fabs(task.xEnd)) / task.nX;

    task.timeStepY = (fabs(task.xStart) + fabs(task.xEnd)) / task.nY;
    task.timeStepZ = (fabs(task.xStart) + fabs(task.xEnd)) / task.nZ;
}