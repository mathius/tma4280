// TMA4280 Supercomputing, Introduction
// Martin Ukrop (martiu@stud.ntnu.no)
// problem set 04
// 2014-02-11

#include "ps04.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <cmath>

// macro for minimum
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

int main(int argc, char** argv) {
    initialize(argc, argv);

    const double pi = 4 * atan(1);
    const double preciseSum = pi*pi/6;

    for (unsigned int i = (int) pow(2,3); i < pow(2,15); i *= 2) {
        double computedSum = computeSum(i);
        if (worldRank == 0) {
            std::cout << "n: " << std::setw(7) << std::left << i << "diff: " << std::setprecision(20) << (preciseSum - computedSum) << std::endl;
        }
    }

    finalize();
    return EXIT_SUCCESS;
}

double computeSum(int vectorLength) {
    double* generatedVector = NULL;
    double* receivedData = NULL;

    // generate vector
    // only done by process 0
    if (worldRank == 0) {
        generatedVector = new double[vectorLength];
#pragma omp parallel for schedule(static)
        for (int i = 0; i < vectorLength; i++) {
            generatedVector[i] = (double) 1 / ((i+1)*(i+1));
        }
    }

// partition vector among processes
#ifdef HAVE_MPI
    receivedData = new double[getLocalLength(vectorLength,worldRank)];
    int* counts = new int[worldSize];
    int* displacements = new int[worldSize];
    // generate lengths and displacements
#pragma omp parallel for schedule(static)
    for (int i = 0; i < worldSize; i++) {
        counts[i] = getLocalLength(vectorLength, i);
        displacements[i] = i * (vectorLength / worldSize) + MIN(vectorLength % worldSize, i);
    }
    MPI::COMM_WORLD.Scatterv(generatedVector, counts, displacements, MPI::DOUBLE, receivedData, counts[worldRank], MPI::DOUBLE, 0);
    delete[] counts;
    delete[] displacements;
#else
    receivedData = generatedVector;
    generatedVector = NULL;
#endif
    delete[] generatedVector;

    // compute local sum (starting from smaller numbers)
    double sum = 0;
#pragma omp parallel for schedule(static) reduction(+:sum)
    for (int i = getLocalLength(vectorLength,worldRank)-1; i >= 0; i--) {
        sum += receivedData[i];
    }

    // reduce if needed (= if using MPI)
#ifdef HAVE_MPI
    double tmpSum = sum;
    MPI::COMM_WORLD.Reduce(&tmpSum, &sum, 1, MPI::DOUBLE, MPI::SUM, 0);
#endif

    delete[] receivedData;
    return sum;
}

int getLocalLength(int globalLength, int rank) {
    int length = globalLength / worldSize;
    if (rank < globalLength % worldSize) {
        length++;
    }
    return length;
}

void initialize(int argc, char** argv) {
#ifdef HAVE_MPI
    MPI::Init(argc, argv);
    worldSize = MPI::COMM_WORLD.Get_size();
    worldRank = MPI::COMM_WORLD.Get_rank();
#else
    worldSize = 1;
    worldRank = 0;
#endif
}

void finalize() {
#ifdef HAVE_MPI
    MPI::Finalize();
#endif
}
