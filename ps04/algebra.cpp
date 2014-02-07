#include "algebra.h"
#include <sys/time.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#ifdef HAVE_OPENMP
#include "omp.h"
#endif

Vector::Vector(int length, double* existingData) {
    lengthGlobal = length;
    commSize = getWorldSize();
    commRank = getWorldRank();
    lengthLocal = lengthGlobal / commSize;
    if (commRank < lengthGlobal % commSize) {
        lengthLocal++;
    }
    offsetLocal = commRank * (lengthGlobal / commSize) + fmin(lengthGlobal % commSize, commRank);
    if (existingData == NULL) {
        data = new double[lengthLocal];
    } else {
        data = existingData;
    }
}

Vector::~Vector() {
    delete[] data;
}

void initialize(int argc, char** argv) {
#ifdef HAVE_MPI
    MPI::Init(argc, argv);
#endif
}

void finalize() {
#ifdef HAVE_MPI
    MPI::Finalize();
#endif
}

int getWorldRank() {
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

double wallTime() {
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
