//
// Created by lenferd on 09.03.17.
//

#include "Task.h"

int Task::initTaskUsingFile(string filename) {
    FILE *infile = fopen(filename.c_str(), "r");

    if (infile == NULL) {
        printf("File reading error. Try to relocate input file\n");
        exit(0);
    }

    // File reading
    fscanf(infile, "XSTART=%lf\n", &xStart);    // start coordinate
    fscanf(infile, "XEND=%lf\n", &xEnd);        // end coordinate
    fscanf(infile, "SIGMA=%lf\n", &sigma);      // coef of heat conduction
    fscanf(infile, "NX=%d\n", &nX);             // count of initial elements?
    fscanf(infile, "TSTART=%lf\n", &tStart);    // start time
    fscanf(infile, "TFINISH=%lf\n", &tFinish);   // finish time
    fscanf(infile, "dt=%lf\n", &dt);            // delta of time difference
    fscanf(infile, "BC=%d\n", &bc);         // Not using right now

    initMemory();

    // Read file
    for (int i = 1; i <= nX; i++) {
        fscanf(infile, "%lf\n", &vect[0][i]);
    }
    fclose(infile);
    return 0;
}

Task::Task() {}

void Task::printTaskData() {
    printf("xStart %lf; xEnd %lf; sigma %lf; nX %d; tStart %lf; tFinish %lf; dt %.10lf;\n",
           xStart, xEnd, sigma, nX, tStart, tFinish, dt);
}

int Task::initMemory(int size) {
    double** vect = new double*[2]; // Vect prev and now
    vect[0] = new double[size + 2]; // +2 because we don't init boundary condition
    vect[1] = new double[size + 2];

    return 0;
}
