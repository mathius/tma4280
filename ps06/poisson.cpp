/*
 * TMA4280 Supercomputing, Introduction
 * Problem set 6
 *
 * C-program to solve the two-dimensional Poisson equation on
 * a unit square using one-dimensional eigenvalue decompositions
 * and fast sine transforms
 *
 * einar m. ronquist
 * ntnu, october 2000
 * revised, october 2001
 *
 * paralelized by Martin Ukrop, 2014
*/

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <memory.h>
#include <cmath>
#include <sys/time.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#ifdef HAVE_OPENMP
#include "omp.h"
#endif

using namespace std;

typedef double Real;

int worldSize;
int worldRank;
int threadMax;
int threadNum;

typedef struct Matrix {
    int dimension;
    int* numColumns;
    int* colDisplacements;
    Real* rawData;
    Real** cols;
    Matrix(int n) {
        dimension = n;
        numColumns = new int[worldSize];
        for (int i = 0; i < worldSize; i++) {
            numColumns[i] = dimension / worldSize;
            if (worldSize-worldRank < dimension%worldSize) {
                numColumns[i]++;
            }
        }
        colDisplacements = new int[worldSize];
        colDisplacements[0] = 0;
        for (int i = 1; i < worldSize; i++) {
            colDisplacements[i] = colDisplacements[i-1]+numColumns[i-1];
        }
        rawData = new Real[dimension*numColumns[worldRank]];
        cols = new Real*[numColumns[worldRank]];
        for (int i = 0; i < worldSize; i++) {
            cols[i] = &(rawData[i*dimension]);
        }
    }
    ~Matrix() {
        delete[] cols;
        delete[] rawData;
        delete[] colDisplacements;
        delete[] numColumns;
    }
} Matrix;

// Sine Transform function prototypes (Fortran, use calling as C)
extern "C" {
    void fst_(Real *v, int *n, Real *w, int *nn);
    void fstinv_(Real *v, int *n, Real *w, int *nn);
}

Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
void transpose (Real **bt, Real **b, int m);
bool parseArguments(int argc, char** argv, int& n);
void initialize(int argc, char** argv);
void finalize();
double wallTime();

// the total number of grid points in each spatial direction is (n+1)
// the total number of degrees-of-freedom in each spatial direction is (n-1)
// this version requires n to be a power of 2

Real exact(Real x, Real y) {
  return x*(pow(x,5)-1.0)*y*(pow(y,5)-1.0);
}

Real source(Real x, Real y) {
  return -30.0*pow(y,4)*x*(pow(x,5.0)-1)-30.0*pow(x,4)*y*(pow(y,5)-1);
}

int main(int argc, char **argv ) {
    // initialize MPI
    initialize(argc, argv);

    Real *diag, **b, **bt, **z, **exactSolution;
    Real pi, h, umax;
    int i, j, n, m, nn;

    if (!parseArguments(argc, argv, n)) {
        return EXIT_FAILURE;
    }

    m  = n-1;
    nn = 4*n;

    diag = createRealArray (m);
    b    = createReal2DArray (m,m);
    bt   = createReal2DArray (m,m);
    z    = createReal2DArray (threadMax,nn);
    exactSolution = createReal2DArray (m,m);

    h    = 1./(Real)n;
    pi   = 4.*atan(1.);

#pragma omp parallel for schedule(static) private(i)
    for (j=0; j < m; j++) {
        diag[j] = 2.*(1.-cos((j+1)*pi/(Real)n));
        for (i=0; i < m; i++) {
            b[j][i] = h*h*source((j+1.0)/n,(i+1.0)/n);
            exactSolution[j][i] = exact((j+1.0)/n,(i+1.0)/n);
        }
    }

    double startTime = wallTime();

#pragma omp parallel private(threadNum)
    {
#ifdef HAVE_OPENMP
        threadNum = omp_get_thread_num();
#endif
#pragma omp for schedule(static)
        for (j=0; j < m; j++) {
            fst_(b[j], &n, z[threadNum], &nn);
        }

#pragma omp single
        transpose(bt,b,m);

#pragma omp for schedule(static)
        for (i=0; i < m; i++) {
            fstinv_(bt[i], &n, z[threadNum], &nn);
        }

#pragma omp for schedule(static) private(i)
        for (j=0; j < m; j++) {
            for (i=0; i < m; i++) {
                bt[j][i] = bt[j][i]/(diag[i]+diag[j]);
            }
        }

#pragma omp for schedule(static)
        for (i=0; i < m; i++) {
            fst_(bt[i], &n, z[threadNum], &nn);
        }

#pragma omp single
        transpose(b,bt,m);

#pragma omp for schedule(static)
        for (j=0; j < m; j++) {
            fstinv_(b[j], &n, z[threadNum], &nn);
        }
    }

    double endTime = wallTime();

    umax = 0.0;
    Real tmp;
    for (j=0; j < m; j++) {
        for (i=0; i < m; i++) {
            tmp = fabs(b[j][i]-exactSolution[j][i]);
            if (tmp > umax) umax = tmp;
        }
    }
    cout << "umax = " << umax << endl;
    cout << "time = " << (endTime-startTime) << endl;

    // finalize MPI
    finalize();
    return EXIT_SUCCESS;
}

void transpose (Real **bt, Real **b, int m)
{
    int i, j;
    for (j=0; j < m; j++) {
        for (i=0; i < m; i++) {
            bt[j][i] = b[i][j];
        }
    }
}

Real *createRealArray (int n)
{
    Real *a;
    int i;
    a = (Real *)malloc(n*sizeof(Real));
    for (i=0; i < n; i++) {
        a[i] = 0.0;
    }
    return (a);
}

Real **createReal2DArray (int n1, int n2)
{
    int i, n;
    Real **a;
    a    = (Real **)malloc(n1   *sizeof(Real *));
    a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
    for (i=1; i < n1; i++) {
        a[i] = a[i-1] + n2;
    }
    n = n1*n2;
    memset(a[0],0,n*sizeof(Real));
    return (a);
}

/**
 * @brief check the second argument is a positive power of 2 and parse it to n
 * @param argc  num of arguments
 * @param argv  cli arguments
 * @param n     parsed problem size
 * @return      parsing process success
 */
bool parseArguments(int argc, char** argv, int& n) {
    if (argc < 2) {
        cerr << "error: Need a problem size." << endl;
        return false;
    }
    istringstream ss(argv[1]);
    if (!(ss >> n)) {
        cerr << "error: Could not parse problem size from '" << string(argv[1]) << "'." << endl;
        return false;
    }
    if (n <= 0 || ( n&(n-1) != 0)) {
        cerr << "error: Problem size of " << n << " is unacceptable (need positive power of 2)." << endl;
        return false;
    }
    return true;
}

/**
 * @brief initialize MPI/OpenMP, fill global variables, print debug info
 * @param argc  num of arguments
 * @param argv  cli arguments
 */
void initialize(int argc, char** argv) {
#ifdef HAVE_MPI
    MPI::Init(argc, argv);
    worldSize = MPI::COMM_WORLD.Get_size();
    worldRank = MPI::COMM_WORLD.Get_rank();
#else
    worldSize = 1;
    worldRank = 0;
#endif
#ifdef HAVE_OPENMP
    threadMax = omp_get_max_threads();
#else
    numThreads = 1;
#endif
    threadNum = 0;
    if (worldRank == 0) {
        cout << "nodes: " << worldSize << endl;
        cout << "max-threads: " << threadMax << endl;
    }
}

/**
 * @brief finalize MPI
 */
void finalize() {
#ifdef HAVE_MPI
    MPI::Finalize();
#endif
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
