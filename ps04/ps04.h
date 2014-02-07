#ifndef PS04_H
#define PS04_H

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#ifdef HAVE_OPENMP
#include "omp.h"
#endif

int main(int argc, char** argv);
double computeSum(int vectorLength);

int getLocalLength(int globalLength, int rank);
int getLocalOffset(int globalLength, int rank);

void initialize(int argc, char** argv);
void finalize();
int getMyWorldRank();
int getWorldSize();

//! \brief Get current wall-clock time
//! \return The current wall time in seconds
double WallTime();

#endif
