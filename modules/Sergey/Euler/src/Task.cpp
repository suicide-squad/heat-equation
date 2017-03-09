//
// Created by lenferd on 09.03.17.
//

#include "Task.h"

Task::Task() {}

int Task::initTaskUsingFile(string filename) {
    FILE *infile = fopen(filename.c_str(), "r");

    if (infile == NULL) {
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
    fscanf(infile, "XSTART=%lf\n", &xStart);    // start coordinate
    fscanf(infile, "XEND=%lf\n", &xEnd);        // end coordinate

    fscanf(infile, "YSTART=%lf\n", &yStart);    // start coordinate
    fscanf(infile, "YEND=%lf\n", &yEnd);        // end coordinate

    fscanf(infile, "ZSTART=%lf\n", &zStart);    // start coordinate
    fscanf(infile, "ZEND=%lf\n", &zEnd);        // end coordinate

    fscanf(infile, "SIGMA=%lf\n", &sigma);      // coef of heat conduction

    fscanf(infile, "NX=%d\n", &nX);             // count of initial elements
    fscanf(infile, "NY=%d\n", &nY);             //
    fscanf(infile, "NZ=%d\n", &nZ);             //

    fscanf(infile, "TSTART=%lf\n", &tStart);    // start time
    fscanf(infile, "TFINISH=%lf\n", &tFinish);   // finish time
    fscanf(infile, "dt=%lf\n", &dt);            // delta of time difference
    fscanf(infile, "BC=%d\n", &bc);         // Not using right now

    initMemory();

    // TODO Нужно писать обращение к элементу, чтобы нормально тут считать файл. Перегрузка?
    // Read file
  /*  for (int i = 1; i <= nX; i++) {
        fscanf(infile, "%lf\n", &vect[0][i]);
    }
    fclose(infile);*/


    return 0;
}

void Task::printTaskData() {
    printf("xStart %lf; xEnd %lf; sigma %lf; nX %d; tStart %lf; tFinish %lf; dt %.10lf;\n",
           xStart, xEnd, sigma, nX, tStart, tFinish, dt);
}

int Task::initMemory() {
    double** vect = new double*[2]; // Vect prev and now

    // TODO Проверить
    int size = (nX + 2) * (nY + 2) * (nZ + 2); // +2 because we don't init boundary condition

    vect[0] = new double[size];
    vect[1] = new double[size];

    return 0;
}

void Task::preparationData() {
    stepSize = (fabs(xStart) + fabs(xEnd)) / nX;

    prevTime = 0;
    currTime = 1;

    vect[0][0] = vect[0][1];
    vect[0][nX+1] = vect[0][nX];

//    double timesize = (tFinal - tStart) / dt;

}
