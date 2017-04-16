#include <iostream>
#include <mpi.h>
#include "Task.h"
#include "SparseMatrix.h"

using std::string;

double getVectorValue(double *vect, int x, int y, int z, Task task);
int main(int argc, char** argv) {

    // Timing variables
    double startTime, endTime;
    const int ROOT = 0;
//    const int ADD_CELL = 2;  // Additional cell (in each size)

    int sizeP = 0, rankP = 0;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankP);

    Task task;

    double *vect;

    if (rankP == ROOT) {
        // File variables
        string functionFile = "../../initial/function.txt";
        string settingFile = "../../initial/setting.ini";

        // Read task settings
        initTaskUsingFile(task, settingFile);
        setTimestep(task);

        // Init memory & read function file
        initMemoryReadDataMPI(vect, functionFile, task);
    }

    MPI_Bcast(&task.xStart, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.xEnd, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.yStart, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.yEnd, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.zStart, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.zEnd, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.sigma, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.bc, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.timeStepX, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.timeStepY, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.timeStepZ, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.nX, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.nY, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.nZ, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.fullVectSize, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.tStart, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.tFinish, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.dt, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);



    int realSizeX = task.nX + 2;
    int realSizeY = realSizeX;
    int realSizeZ = realSizeY * task.nY;

    int lineSizeP = sizeP / 2;
    int proc_nX = realSizeX;
    int proc_nY = task.nY / lineSizeP;
    int proc_nZ = task.nZ / lineSizeP;
    int proc_vect_size = proc_nX * proc_nY * proc_nZ;
    int block_size = proc_nX * proc_nY;
//    int block_size_add = block_size + ADD_CELL * 2;

    // Generate data for send
    double *sendvect = new double [sizeP * proc_vect_size];
    for (int l = 0; l < sizeP; ++l) {
        for (int k = 0; k < proc_nZ; ++k) {
            for (int i = 0; i < block_size; ++i) {
                sendvect[l*proc_vect_size + k * block_size + i]
                        = vect[proc_nZ * realSizeZ * l + k * block_size+ i];
            }
        }
    }


    int *displs = new int[sizeP];
    int *sendcounts = new int[sizeP];


    // MAY NOT WORK WELL
    for (int l = 0; l < sizeP; ++l) {
        displs[l] = proc_nZ * realSizeZ * l;
        sendcounts[l] = realSizeY * proc_nY;     // we send half part
    }

    // Sending should be in two steps

    double *proc_vect = new double[proc_vect_size];

    // STEP 1
    MPI_Scatterv(vect, sendcounts, displs, MPI_DOUBLE,
                 proc_vect, proc_vect_size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);



    // Prepare STEP 2

    // MAY NOT WORK WELL
    for (int l = 0; l < sizeP; ++l) {
        displs[l] = proc_nZ * realSizeZ * l;
        sendcounts[l] = realSizeY * proc_nY;     // we send half part
    }

    // STEP 1
    MPI_Scatterv(vect, sendcounts, displs, MPI_DOUBLE,
                 proc_vect, proc_vect_size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);




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
    FILE *outfile = fopen("../../result/Sergey/Sergey_Euler.txt", "w");

    double outData;
    for (int i = 1; i <= task.nX; ++i) {
        fprintf(outfile, "%2.15le\n", getVectorValue(vect[0],i,0,0,task));

    }
}

double getVectorValue(double *vect, int x, int y, int z, Task task) {
    return vect[x + (task.nX + 2) * y + (task.nX+2)*task.nY*z];
}

