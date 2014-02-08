// TMA4280 Supercomputing, Introduction
// Martin Ukrop (martiu@stud.ntnu.no)
// problem set 04
// 2014-02-08

#ifndef PS04_H
#define PS04_H

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#ifdef HAVE_OPENMP
#include "omp.h"
#endif

// global variables, more efficient than asking MPI every time
int worldSize;
int worldRank;

// main program, runs computeSum in a loop, produces output
int main(int argc, char** argv);

// computes S_n = Sum_{i=1}^{vectorLength} (1/(i^2))
double computeSum(int vectorLength);

// returns part of globalLength partitioned to process rank in World
int getLocalLength(int globalLength, int rank);

// initialize MPI, set worldSize and worldRank
void initialize(int argc, char** argv);

// finalize MPI
void finalize();

#endif
