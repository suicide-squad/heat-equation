#include <iostream>
#include <mpi.h>
#include "Task.h"
#include "SparseMatrix.h"

using std::string;
using std::cout;
using std::endl;
using std::flush;

double getVectorValue(double *vect, int x, int y, int z, Task task);
int main(int argc, char** argv) {

    // Timing variables
    double startTime, endTime;
    const int ROOT = 0;
//    const int ADD_CELL = 2;  // Additional cell (in each size)

    int sizeP, rankP;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankP);

    Task task;

    double *vect;

    if (rankP == ROOT) {
        // File variables
        string functionFile = "../../initial/function.txt";
//        string functionFile = "../../initial_test/function.txt";
        string settingFile = "../../initial/setting.ini";
//        string settingFile = "../../initial_test/setting.ini";

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


    if (sizeP % 4 != 0) {
        cout << "Error! Wrong proc count" << endl;
        exit(0);
    }

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

    int proc_realSizeY = realSizeX;
    int proc_realSizeZ = proc_realSizeY * proc_nY;


    // Generate data for send
    double *sendvect;

    if (rankP == ROOT) {
        sendvect = new double [sizeP * proc_vect_size];

        for (int i = 0; i < sizeP * proc_vect_size; ++i) {
            sendvect[i] = -1000;
        }

        cout << "sizeP: " << sizeP << endl;
        cout << "proc_nX: " << proc_nX << endl;
        cout << "proc_nY: " << proc_nY << endl;
        cout << "proc_nZ: " << proc_nZ << endl;
        cout << "block_size: " << block_size << endl;
        cout << "proc_vect_size: " << proc_vect_size << endl;

//        cout << ": " <<  << endl;

//        startTime = MPI_Wtime();

        for (int p = 0; p < sizeP; ++p) {
            int proc_lvl_offset = (p / lineSizeP) * lineSizeP * proc_nZ * proc_realSizeZ;
            int fist_second_offset = (p % lineSizeP) * proc_nY * proc_realSizeY;
            for (int z = 0; z < proc_nZ; ++z) {
                for (int i = 0; i < block_size; ++i) {
                    sendvect[p * proc_vect_size + z * block_size + i] =
                            vect[proc_lvl_offset + fist_second_offset + z * realSizeZ + i];
                }
            }
        }

//        endTime = MPI_Wtime();
//        printf("Repack run time %.15lf\n", endTime - startTime);
    }

//    if (rankP == ROOT) {
//        for (int i = 0; i < sizeP * proc_vect_size; ++i) {
//            cout << i << "  " << sendvect[i] << endl;
//        }
//    }


    int *displs = new int[sizeP];
    int *sendcounts = new int[sizeP];

    double **proc_vect = new double*[2];
    proc_vect[0] = new double[proc_vect_size];
    proc_vect[1] = new double[proc_vect_size];

    for (int l = 0; l < sizeP; ++l) {
        displs[l] = proc_vect_size * l;
        sendcounts[l] = proc_vect_size;     // we send half part
    }


    MPI_Scatterv(sendvect, sendcounts, displs, MPI_DOUBLE,
                 proc_vect[0], proc_vect_size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);


//    if (rankP == 3) {
//        cout << "print vectors ========" << endl;
//        for (int i = 0; i < proc_vect_size; ++i) {
//            cout << proc_vect[0][i] << endl << flush;
//        }
//    }

    MPI_Barrier(MPI_COMM_WORLD);
    // value for the matrix
    MatrixValue matrixValue;
    matrixValue.x1 = (task.sigma * task.dt) / (task.timeStepX * task.timeStepX);
    matrixValue.y1 = (task.sigma * task.dt) / (task.timeStepY * task.timeStepY);
    matrixValue.z1 = (task.sigma * task.dt) / (task.timeStepZ * task.timeStepZ);
    matrixValue.x2Comp = (1 - 2 * matrixValue.x1 - 2 * matrixValue.y1 - 2 * matrixValue.z1);

//    cout << "matrixValue.x1: " << matrixValue.x1 << endl;
//    cout << "matrixValue.y1: " << matrixValue.y1 << endl;
//    cout << "matrixValue.z1: " << matrixValue.z1 << endl;
//    cout << "matrixValue.x2C: " << matrixValue.x2Comp << endl;

    // init and fill sparseMatrix
    SparseMatrix spMat;
    int sparseMatrixSize = 9 * proc_nX * proc_nY * proc_nZ;

    spMatrixInit(spMat, sparseMatrixSize, proc_vect_size);

//    cout << "proc_vect_size: " << proc_vect_size << endl;
    // КОСТЫЛИИИИИИИИИИИИИИИИИи
    fillMatrix3d6Expr_wo_boundaries(spMat, matrixValue, proc_nX - 2, proc_nY, proc_nZ);


    proc_realSizeY = realSizeX;
    proc_realSizeZ = proc_realSizeY * proc_nY;

    int prevTime = 0;
    int currTime = 1;

    startTime = MPI_Wtime();

    int send_vect_size = proc_nZ * proc_nX; // Because nY don't involved
    double *left_vect = new double[send_vect_size];
    double *right_vect = new double[send_vect_size];

    for (double j = 0; j < task.tFinish; j += task.dt) {
        multiplicateVector(spMat, proc_vect[prevTime], proc_vect[currTime], proc_vect_size);

        // After each iteration swap data

        // 1 and 2 iter - IT'S TIME TO REPACK
        for (int z = 0; z < proc_nZ; ++z) {
            for (int i = 0; i < proc_nX; ++i) {
                left_vect[z * proc_nX + i] = proc_vect[currTime][z * proc_realSizeZ + i];
                right_vect[z * proc_nX + i]
                        = proc_vect[currTime][z * proc_realSizeZ + (proc_nY - 1) * proc_realSizeY + i];
            }
        }

        int left_proc = (rankP / lineSizeP) * lineSizeP + (rankP - 1 + lineSizeP) % lineSizeP;
        int right_proc = (rankP / lineSizeP) * lineSizeP + (rankP + 1 + lineSizeP) % lineSizeP;

        MPI_Sendrecv(left_vect, send_vect_size, MPI_DOUBLE, left_proc, 0,
                     left_vect, send_vect_size, MPI_DOUBLE, right_proc, 0, MPI_COMM_WORLD, &status);

        MPI_Sendrecv(right_vect, send_vect_size, MPI_DOUBLE, right_proc, 0,
                     right_vect, send_vect_size, MPI_DOUBLE, left_proc, 0, MPI_COMM_WORLD, &status);

        for (int z = 0; z < proc_nZ; ++z) {
            for (int i = 0; i < proc_nX; ++i) {
                proc_vect[currTime][z * proc_realSizeZ + i] = left_vect[z * proc_nX + i];
                proc_vect[currTime][z * proc_realSizeZ + (proc_nY - 1) * proc_realSizeY + i] = right_vect[z * proc_nX + i];
            }
        }

        int top_proc = (rankP - lineSizeP + sizeP) % sizeP;
        int bottom_proc = (rankP + lineSizeP) % sizeP;

//        cout << "top_proc: " << top_proc << endl;
//        cout << "bottom_proc: " << bottom_proc << endl;

        // 1 and 2 iter - in one field in the memory
        send_vect_size = proc_nY * proc_nX;
        // top
        MPI_Sendrecv(proc_vect[currTime], send_vect_size, MPI_DOUBLE, top_proc, 0,
                     proc_vect[currTime], send_vect_size, MPI_DOUBLE, bottom_proc, 0, MPI_COMM_WORLD, &status);

        // bottom
        int offset = proc_vect_size - send_vect_size;
        MPI_Sendrecv(proc_vect[currTime] + offset, send_vect_size, MPI_DOUBLE, bottom_proc, 0,
                     proc_vect[currTime] + offset, send_vect_size, MPI_DOUBLE, top_proc, 0, MPI_COMM_WORLD, &status);


        int boundaries_offset = proc_realSizeY;
        // BOUNDARIES!
        for (int k = 0; k < proc_nY * proc_nZ; ++k) {
            proc_vect[currTime][k * boundaries_offset] = proc_vect[currTime][k * boundaries_offset + 1];
            proc_vect[currTime][k * boundaries_offset + proc_nX - 1] =
                    proc_vect[currTime][k * boundaries_offset + proc_nX - 2];
        }
        prevTime = (prevTime + 1) % 2;
        currTime = (currTime + 1) % 2;
    }

    // And now we should scatter it in one vector...

    MPI_Gather(proc_vect[0], proc_vect_size, MPI_DOUBLE,
               sendvect, proc_vect_size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

//    if (rankP == ROOT) {
//        cout << "print sendvect vect" << endl;
//        cout.precision(15);
//        for (int i = 0; i < sizeP * proc_vect_size; ++i) {
//            cout << i << "  " << sendvect[i] << endl;
//        }
//        cout << "end" << endl;
//    }in_Vect_offset



    if (rankP == ROOT) {
        for (int j = 0; j < task.fullVectSize; ++j) {
            vect[j] = -10000;
        }

        cout << "sizeP: " << sizeP << endl;
        cout << "proc_nX: " << proc_nX << endl;
        cout << "proc_nY: " << proc_nY << endl;
        cout << "proc_nZ: " << proc_nZ << endl;
        cout << "block_size: " << block_size << endl;
        cout << "proc_vect_size: " << proc_vect_size << endl;
        cout << "realSizeZ: " << realSizeZ << endl;
        cout << "lineSizeP: " << lineSizeP << endl;
//        cout << ": " <<  << endl;

        for (int p = 0; p < sizeP; ++p) {
            int proc_lvl_offset = (p / lineSizeP) * lineSizeP * proc_nZ * proc_realSizeZ;
            int fist_second_offset = (p % lineSizeP) * proc_nY * proc_realSizeY;
            for (int z = 0; z < proc_nZ; ++z) {
                for (int i = 0; i < block_size; ++i) {
                    vect[proc_lvl_offset + fist_second_offset + z * realSizeZ + i] =
                            sendvect[p * proc_vect_size + z * block_size + i];
                }
            }
        }
    }


//    for (int l = 0; l < sizeP; ++l) {
//        for (int k = 0; k < proc_nZ; ++k) {
//            for (int i = 0; i < block_size; ++i) {
//                sendvect[l * proc_vect_size + k * block_size + i]
//                        = vect[(l % lineSizeP) * lineSizeP * realSizeZ + k * realSizeZ + i];


    if (rankP == ROOT) {
        endTime = MPI_Wtime();
        printf("Run time %.15lf\n", endTime - startTime);


        // Output
        FILE *outfile = fopen("../../result/Sergey/Sergey_Euler_MPI_full.txt", "w");
//        FILE *outfile = fopen("../../result/Sergey/Sergey_Euler_MPI_full_test2.txt", "w");

        for (int i = 0; i < task.fullVectSize; ++i) {
            if ( i % (task.nX + 2) != 0 && i % (task.nX + 2) != task.nX + 1)
                fprintf(outfile, "%2.15le\n", vect[i]);
        }
//        for (int i = 1; i <= task.nX; ++i) {
//            fprintf(outfile, "%2.15le\n", getVectorValue(vect, i, 0, 0, task));
//
//        }
    }
    MPI_Finalize();
}

double getVectorValue(double *vect, int x, int y, int z, Task task) {
    return vect[x + (task.nX + 2) * y + (task.nX+2)*task.nY*z];
}

