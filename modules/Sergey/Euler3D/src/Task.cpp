//
// Created by lenferd on 09.03.17.
//

#include "Task.h"

Task::Task() {}

int Task::initTaskUsingFile(string settingFile, string functionFile) {
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
    fscanf(inSettingfile, "XSTART=%lf\n", &xStart);    // start coordinate
    fscanf(inSettingfile, "XEND=%lf\n", &xEnd);        // end coordinate

    fscanf(inSettingfile, "YSTART=%lf\n", &yStart);    // start coordinate
    fscanf(inSettingfile, "YEND=%lf\n", &yEnd);        // end coordinate

    fscanf(inSettingfile, "ZSTART=%lf\n", &zStart);    // start coordinate
    fscanf(inSettingfile, "ZEND=%lf\n", &zEnd);        // end coordinate

    fscanf(inSettingfile, "SIGMA=%lf\n", &sigma);      // coef of heat conduction

    fscanf(inSettingfile, "NX=%d\n", &nX);             // count of initial elements
    fscanf(inSettingfile, "NY=%d\n", &nY);             //
    fscanf(inSettingfile, "NZ=%d\n", &nZ);             //

    fscanf(inSettingfile, "TSTART=%lf\n", &tStart);    // start time
    fscanf(inSettingfile, "TFINISH=%lf\n", &tFinish);   // finish time
    fscanf(inSettingfile, "dt=%lf\n", &dt);            // delta of time difference
    fscanf(inSettingfile, "BC=%d\n", &bc);         // Not using right now

    fclose(inSettingfile);
    FILE *inFunctionfile = fopen(functionFile.c_str(), "r");

    initMemory();

    double datakek;
    /// Read file
    for (int k = 0; k < nZ; k++) {
        for (int j = 0; j < nY; ++j) {
            for (int i = 1; i < nX + 1; ++i) {
                fscanf(inFunctionfile, "%lf\n", &vect[0][i + (nX + 2) * j + (nX+2)*nY*k]);
            }
        }
    }

    fclose(inFunctionfile);
    return 0;
}

void Task::printTaskData() {
    printf("xStart %lf; xEnd %lf; sigma %lf; nX %d; tStart %lf; tFinish %lf; dt %.10lf;\n",
           xStart, xEnd, sigma, nX, tStart, tFinish, dt);
}

int Task::initMemory() {
    vect = new double*[2]; // Vect prev and now

    fullVectSize = (nX + 2) * (nY) * (nZ);  // +2 because we don't init boundary condition
                                        // And that only on X coordinate

    vect[0] = new double[fullVectSize];
    vect[1] = new double[fullVectSize];

    for (int i = 0; i < fullVectSize; ++i) {
        vect[0][i] = 0;
        vect[1][i] = 0;
    }
    return 0;
}

void Task::preparationData() {
    timeStepX = (fabs(xStart) + fabs(xEnd)) / nX;

    timeStepY = (fabs(xStart) + fabs(xEnd)) / nY;
    timeStepZ = (fabs(xStart) + fabs(xEnd)) / nZ;

    prevTime = 0;
    currTime = 1;

    vect[0][0] = vect[0][1];
    vect[0][nX+1] = vect[0][nX];

//    double timesize = (tFinal - tStart) / dt;

}

void Task::setData(double data, int x, int y, int z, int timeVect) {
    vect[timeVect][x + (nX + 2) * y + (nX+2)*nY*z] = data;
}

double Task::getData(int x, int y, int z, int timeVect) {
    return vect[timeVect][x + (nX + 2) * y + (nX+2)*nY*z];
}

void Task::printVect(int timeVect) {
    for (int i = 0; i < fullVectSize; ++i) {
        printf("%lf", vect[timeVect][i]);
    }
}

void Task::printVectFile(string filename, int timeVect) {
    FILE *outfile = fopen(filename.c_str(), "w");

    for (int i = 0; i < fullVectSize; ++i) {
        fprintf(outfile, "%lf\n", vect[timeVect][i]);
    }
}

void Task::createMatrix(TaskExpressions &taskexpr) {
    spMatrixInit(matrix, 7*nX*nY*nZ + 2*nY*nZ, fullVectSize);
    fillMatrix3d6Expr(matrix, taskexpr, nX, nY, nZ);
}

void Task::setTaskExpr(TaskExpressions &task) {
    task.x1 = (sigma * dt) / (timeStepX * timeStepX);
    printf("sigma %lf\n", sigma);
    printf("dt %lf\n", dt);
    printf("timeStepX %lf\n", timeStepX);
    printf("taskx1 %lf\n", task.x1);
    task.y1 = (sigma * dt) / (timeStepY * timeStepY);
    printf("tasky1 %lf\n", task.y1);
    task.z1 = (sigma * dt) / (timeStepZ * timeStepZ);
    printf("taskz1 %lf\n", task.z1);

    task.x2Comp = (1 - 2 * task.x1 - 2 * task.y1 - 2 * task.z1);
    printf("taskx2 %lf\n", task.x2Comp);

    //1 - 2 * ((sigma * dt) * 1 / (timeStepX * timeStepX) * 1 / (timeStepY * timeStepY) * 1/ (timeStepZ * timeStepZ);

}

