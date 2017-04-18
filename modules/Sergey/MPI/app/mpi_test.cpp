//
// Created by lenferd on 14.04.17.
//

#include "mpi.h"
#include <istream>

int main(int argc, char **argv) {
    int sizeP;
    int rankP;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankP);

    std::cout << "proc ranK " << rankP << " of " << sizeP << std::endl;

    MPI_Finalize();
    return 0;
}