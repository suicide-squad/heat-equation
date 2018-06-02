//
// Created by kirill on 19.03.17.
//

// Implicit Euler Method. Method Jacobi

#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include <stdbool.h>

#include <fcntl.h>
#include <unistd.h>

#include "utils/ts.h"
#include "multMV.h"
#include "createSpMat.h"
#include "parser.h"
#include "sgpu.h"

#define EPS 1e-14

#define DIM_CART 2
#define SHIFT 1
#define RESERVE SHIFT*2

#define IND(x, y, z) ((x) + (y)*NX + (z)*NX*(NYr + RESERVE))

static inline bool dist(const TYPE * const x1, const TYPE * const x2, const int N) {
    for (int i = 0; i < N; i++)
        if (ABS(x1[i] - x2[i]) > EPS)
            return false;
    return true;
}

int main(int argc, char **argv) {
    double t0 = 0.0, t1 = 0.0;
    const size_t len = 80;
    char nameHost[len];
    gethostname(nameHost, len);
    
    int sizeP = 1;
    int rankP = ROOT;
    
    size_t sizeTime;
    
    int blockYP = 0, blockZP = 0;
    
    // INIT MPI SETTING
    MPI_Status status[4];
    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankP);
    
    if (rankP == ROOT) {
        printf("MPI RUN. %d size processes\n", sizeP);
        printf("name host - %s\n", nameHost);
        fflush(stdout);
    }
    MPI_Comm gridComm;
    
    SpMatrix mat;

    TYPE coeffs[3];
    TYPE *u = NULL, *u_chunk = NULL, *un_chunk = NULL;
    TYPE k;
    int NX, NY, NZ;
    int dim = 0;
    
    if (rankP == ROOT) {
        Setting setting;
        int error;
    
        error = readSetting(INPUT_EULER_SETTING_PATH, &setting);
        
        if (error != OK) return error;
        
        dim = (setting.NX + RESERVE) * (setting.NY + RESERVE) * (setting.NZ + RESERVE);
        u = (TYPE *) calloc((size_t) dim, sizeof(TYPE));
        
        error = readFunction(INPUT_EULER_FUNCTION_PATH, u, setting.NX + RESERVE,
                             setting.NY + RESERVE, setting.NZ + RESERVE, SHIFT);
        
        if (error != OK) return error;
        
        sizeTime = (size_t) ((setting.TFINISH - setting.TSTART) / setting.dt);
        
#if ENABLE_PARALLEL
        printf("PARALLEL VERSION!\n");
        // omp_set_num_threads(1);
#endif
        const TYPE dx = ABS(setting.XSTART - setting.XEND) / setting.NX;
        const TYPE dy = ABS(setting.YSTART - setting.YEND) / setting.NY;
        const TYPE dz = ABS(setting.ZSTART - setting.ZEND) / setting.NZ;
        
        coeffs[0] = setting.dt * setting.SIGMA / (dx * dx);
        coeffs[1] = setting.dt * setting.SIGMA / (dy * dy);
        coeffs[2] = setting.dt * setting.SIGMA / (dz * dz);
        
        k = 1.0 / (1.0 + 2.0 * setting.SIGMA * setting.dt * (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)));
        
        // coeffs[0] = 1; coeffs[1] = 2; coeffs[2] = 3;
        
        NX = setting.NX + RESERVE;
        NY = setting.NY + RESERVE;
        NZ = setting.NZ + RESERVE;
    }
    
    MPI_Bcast(&sizeTime, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(coeffs, 3, MPI_TYPE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_TYPE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&NX, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&NY, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&NZ, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    
    //  Определения числа процессов в каждом измерении
    get_blocks(&blockYP, &blockZP, sizeP);
    
    if (rankP == ROOT) printf("blockY %d blockZ %d\n", blockYP, blockZP);
    
    const int NYr = (NY - RESERVE) / blockYP;
    const int NZr = (NZ - RESERVE) / blockZP;
    
    const int dimChunk = NX * (NYr + RESERVE) * (NZr + RESERVE);
    u_chunk = (TYPE *) aligned_alloc(ALIGNMENT, sizeof(TYPE) * dimChunk);
    un_chunk = (TYPE *) aligned_alloc(ALIGNMENT, sizeof(TYPE) * dimChunk);
    
    int rank_left, rank_right, rank_down, rank_top;
    
    //  размер каждой размерности
    const int dims[] = {blockZP, blockYP};
    //  наличие циклов в каждой размерности
    const int periods[] = {0, 0};
    
    //  разрешение системе менять номера процессов
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, DIM_CART, dims, periods, reorder, &gridComm);
    
    scatter_by_block(u, u_chunk, NX, NY, NYr, NZr, gridComm, RESERVE);
    
    MPI_Cart_shift(gridComm, 1, -1, &rank_left, &rank_right);
    MPI_Cart_shift(gridComm, 0, 1, &rank_down, &rank_top);
    
    const int nonZero = dimChunk * 6u;
    
    if (rankP == ROOT) printf("TimeSize -\t%lu\n", sizeTime);
    
    initSpMat(&mat, nonZero, dimChunk);
    createImplicitSpMat(&mat, coeffs, NX, NYr + RESERVE, NZr + RESERVE);
    
    // Создание типа плоскости XY и XZ
    MPI_Datatype planeXY;
    MPI_Type_vector(NZr + RESERVE, NX, NX * (NYr + RESERVE), MPI_TYPE, &planeXY);
    MPI_Type_commit(&planeXY);
    
    MPI_Datatype planeXZ;
    MPI_Type_contiguous(NX * (NYr + RESERVE), MPI_TYPE, &planeXZ);
    MPI_Type_commit(&planeXZ);
    
//  printSpMat(mat);
    
    TYPE * x1 = (TYPE *) aligned_alloc(ALIGNMENT, sizeof(TYPE) * dimChunk);
    TYPE * x2 = (TYPE *) aligned_alloc(ALIGNMENT, sizeof(TYPE) * dimChunk);
    
    // int countIter = 0;

    copyingBorders(u_chunk, NX, NYr + RESERVE, NZr + RESERVE);


    // ОСНОВНЫЕ ВЫЧИСЛЕНИЯ
    if (rankP == ROOT) {
        printf("START!\n");
        t0 = WTIME();
    }
    for (int t = 1; t <= 100; t++) {

        //    Передача влево по Y
        MPI_Sendrecv(&u_chunk[IND(0, NYr, 0)], 1, planeXY, rank_left, 0,
                     &u_chunk[IND(0, 0, 0)], 1, planeXY, rank_right, 0, gridComm, &status[0]);
        
        //    Передача вправо по Y
        MPI_Sendrecv(&u_chunk[IND(0, 1, 0)], 1, planeXY, rank_right, 1,
                     &u_chunk[IND(0, NYr + 1, 0)], 1, planeXY, rank_left, 1, gridComm, &status[1]);
        
        //    Передача вниз по Z
        MPI_Sendrecv(&u_chunk[IND(0, 0, 1)], 1, planeXZ, rank_down, 2,
                     &u_chunk[IND(0, 0, NZr + 1)], 1, planeXZ, rank_top, 2, gridComm, &status[2]);
        
        //    Передача вверх по Z
        MPI_Sendrecv(&u_chunk[IND(0, 0, NZr)], 1, planeXZ, rank_top, 3,
                     &u_chunk[IND(0, 0, 0)], 1, planeXZ, rank_down, 3, gridComm, &status[3]);
        
        memcpy(x1, u_chunk, dimChunk * sizeof(TYPE));
        
        do {
            
            multMV(x2, mat, x1, NX, (NYr + RESERVE), (NZr + RESERVE), coeffs);
            
            for (int z = 1; z < NZr + RESERVE - 1; z++) {
                for (int y = 1; y < NYr + RESERVE - 1; y++) {
                    for (int x = 1; x < NX - 1; x++) {
                        x2[IND(x, y, z)] = (u_chunk[IND(x, y, z)] + x2[IND(x, y, z)]) * k;
                       if (x == 1)
                           x2[IND(x - 1, y, z)] = x2[IND(x, y, z)];
                       if (x == NX - 2)
                           x2[IND(x + 1, y, z)] = x2[IND(x, y, z)];
                       if (y == 1)
                           x2[IND(x, y - 1, z)] = x2[IND(x, y, z)];
                       if (y == NY - 2)
                           x2[IND(x, y + 1, z)] = x2[IND(x, y, z)];
                       if (z == 1)
                           x2[IND(x, y, z - 1)] = x2[IND(x, y, z)];
                       if (z == NZ - 2)
                           x2[IND(x, y, z + 1)] = x2[IND(x, y, z)];
                    }
                }
            }

            TYPE * tmp = x1;
            x1 = x2;
            x2 = tmp;
            
            // countIter++;
            
        } while (!dist(x1, x2, dimChunk));
        
        // u_chunk = x1;
                memcpy(u_chunk, x1 , dimChunk * sizeof(TYPE));


        
    }
    
    //  *******************
    
    if (rankP == ROOT) {
        printf("FINISH!\n\n");
        t1 = WTIME();
    }
    
    gather_by_block(u, u_chunk, NX, NY, NYr, NZr, RESERVE, gridComm);
    
    if (rankP == ROOT) {
        double diffTime = t1 - t0;
        printf("Time -\t%.3lf\n", diffTime);
        
        writeFunction3D(RESULT_IEULER_PATH, u, NX, NY, NZ, SHIFT);
        
        printf("DONE!!!\n\n");
        free(u);
    }
    
    freeSpMat(&mat);
    free(un_chunk);
    free(x2);
    
    MPI_Type_free(&planeXY);
    MPI_Type_free(&planeXZ);
    MPI_Finalize();
    
    return 0;
}