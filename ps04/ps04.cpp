#include "ps04.h"

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <cmath>

int main(int argc, char** argv) {
    initialize(argc, argv);

    const double pi = 4 * atan(1);
    const double preciseSum = pi*pi/6;

    for (unsigned int i = (int) pow(2,3); i < pow(2,15); i *= 2) {
        double computedSum = computeSum(i);
        if (getMyWorldRank() == 0) {
            std::cout << "n: " << std::setw(7) << std::left << i << "diff: " << preciseSum - computedSum << std::endl;
        }
    }

    finalize();
    return EXIT_SUCCESS;
}

double computeSum(int vectorLength) {

    int localLength = getLocalLength(vectorLength,getMyWorldRank());
    double* receivedData = NULL;

    if (getMyWorldRank() == 0) {
        // generate vector
        double* tmpVector = new double[vectorLength];
#pragma omp parallel for schedule(static)
        for (int i = 0; i < vectorLength; i++) {
            tmpVector[i] = (double) 1 / ((i+1)*(i+1));
        }

#ifdef HAVE_MPI
        // spread vector among processes
#pragma omp parallel for schedule(static)
        for (int i = 1; i < getWorldSize(); i++) {
            MPI::COMM_WORLD.Send(tmpVector + getLocalOffset(vectorLength,i),
                                 getLocalLength(vectorLength,i),
                                 MPI::DOUBLE, i, getLocalLength(vectorLength,i));
        }
#endif
        receivedData = tmpVector;
    }

#ifdef HAVE_MPI
    // recieve vector parts
    if (receivedData == NULL) {
        receivedData = new double[localLength];
        MPI::COMM_WORLD.Recv(receivedData, localLength, MPI::DOUBLE, 0, localLength);
    }
#endif

    // compute sum (starting from smaller numbers)
    double sum = 0;
#pragma omp parallel for schedule(static) reduction(+:sum)
    for (int i = localLength-1; i >= 0; i--) {
        sum += receivedData[i];
    }

    // reduce if needed
#ifdef HAVE_MPI
    double tmpSum = sum;
    MPI::COMM_WORLD.Reduce(&tmpSum, &sum, 1, MPI::DOUBLE, MPI::SUM, 0);
#endif

    delete[] receivedData;
    return sum;
}

int getLocalLength(int globalLength, int rank) {
    int length = globalLength / getWorldSize();
    if (rank < globalLength % getWorldSize()) {
        length++;
    }
    return length;
}

int getLocalOffset(int globalLength, int rank) {
    int offset = rank * (globalLength / getWorldSize());
    offset += fmin(globalLength % getWorldSize(), rank);
    return offset;
}

void initialize(int argc, char** argv) {
#ifdef HAVE_MPI
#ifdef HAVE_OPENMP
  int aquired = MPI::Init_thread(argc, argv, MPI_THREAD_MULTIPLE);
  if (getMyWorldRank() == 0) {
    std::cout << "aquired MPI threading level: ";
    if (aquired == MPI_THREAD_SINGLE)
      std::cout << "MPI_THREAD_SINGLE";
    if (aquired == MPI_THREAD_FUNNELED)
      std::cout << "MPI_THREAD_FUNNELED";
    if (aquired == MPI_THREAD_SERIALIZED)
      std::cout << "MPI_THREAD_SERIALIZED";
    if (aquired == MPI_THREAD_MULTIPLE)
      std::cout << "MPI_THREAD_MULTIPLE";
    std::cout << std::endl;
  }
#else
  MPI::Init(argc, argv);
#endif
#endif
}

void finalize() {
#ifdef HAVE_MPI
    MPI::Finalize();
#endif
}

int getMyWorldRank() {
#ifdef HAVE_MPI
    return MPI::COMM_WORLD.Get_rank();
#endif
    return 0;
}

int getWorldSize() {
#ifdef HAVE_MPI
    return MPI::COMM_WORLD.Get_size();
#endif
    return 1;
}

double WallTime ()
{
#ifdef HAVE_MPI
  return MPI_Wtime();
#elif defined(HAVE_OPENMP)
  return omp_get_wtime();
#else
  struct timeval tmpTime;
  gettimeofday(&tmpTime,NULL);
  return tmpTime.tv_sec + tmpTime.tv_usec/1.0e6;
#endif
}
