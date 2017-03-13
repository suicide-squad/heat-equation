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

    /// Read file
    for (int k = 0; k < nY; k++) {
        for (int j = 0; j < nY; ++j) {
            for (int i = 1; i < nX + 1; ++i) {
                fscanf(infile, "%lf\n", &vect[0][i]);
            }
        }
    }

    fclose(infile);
    return 0;
}

void Task::printTaskData() {
    printf("xStart %lf; xEnd %lf; sigma %lf; nX %d; tStart %lf; tFinish %lf; dt %.10lf;\n",
           xStart, xEnd, sigma, nX, tStart, tFinish, dt);
}

int Task::initMemory() {
    double** vect = new double*[2]; // Vect prev and now

    fullVectSize = (nX + 2) * (nY) * (nZ);  // +2 because we don't init boundary condition
                                        // And that only on X coordinate

    vect[0] = new double[fullVectSize];
    vect[1] = new double[fullVectSize];

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

void Task::setData(double data, int x, int y, int z, int timeVect = 0) {
    vect[timeVect][x + (nY * y) + (nZ * z)] = data;
}

double Task::getData(int x, int y, int z, int timeVect) {
    return vect[timeVect][x + (nY * y) + (nZ * z)];
}

void Task::printVect(int timeVect = 0) {
    for (int i = 0; i < fullVectSize; ++i) {
        printf("%lf", vect[timeVect][i]);
    }
}
