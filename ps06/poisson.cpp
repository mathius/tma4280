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

// the total number of grid points in each spatial direction is (n+1)
// the total number of degrees-of-freedom in each spatial direction is (n-1)
// this version requires n to be a power of 2

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <cmath>
#include <sys/time.h>
#include "poisson.h"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#ifdef HAVE_OPENMP
#include "omp.h"
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::istringstream;

int threadMax;

Matrix::Matrix(int pCols, int pRows, bool mpiDistributed) {
    worldSize = 1;
    worldRank = 0;
#ifdef HAVE_MPI
    if (mpiDistributed) {
        worldSize = MPI::COMM_WORLD.Get_size();
        worldRank = MPI::COMM_WORLD.Get_rank();
    }
#endif
    rows = pRows;
    cols = pCols;
    numCols = new int[worldSize];
    for (int i = 0; i < worldSize; i++) {
        numCols[i] = cols / worldSize;
        // left-over columns are on the last nodes
        if (worldSize-i <= cols%worldSize) {
            numCols[i]++;
        }
    }
    localNumCols = numCols[worldRank];
    sendCounts = new int[worldSize];
    for (int i = 0; i < worldSize; i++) {
        sendCounts[i] = localNumCols*numCols[i];
    }
    colDispl = new int[worldSize];
    colDispl[0] = 0;
    for (int i = 1; i < worldSize; i++) {
        colDispl[i] = colDispl[i-1]+numCols[i-1];
    }
    sendDispl = new int[worldSize];
    sendDispl[0] = 0;
    for (int i = 1; i < worldSize; i++) {
        sendDispl[i] = sendDispl[i-1]+sendCounts[i-1];
    }
    rawData = new double[rows*localNumCols]();
    data = new double*[localNumCols];
    for (int i = 0; i < localNumCols; i++) {
        data[i] = &(rawData[i*rows]);
    }
}
Matrix::~Matrix() {
    delete[] data;
    delete[] rawData;
    delete[] sendCounts;
    delete[] sendDispl;
    delete[] numCols;
    delete[] colDispl;
}

int main(int argc, char **argv ) {
    initialize(argc, argv);
    // problem size (n), general indices (i,j)
    int i, j, n;
    double umax, timeStart, timeEnd;

    if (!parseArguments(argc, argv, n)) { return EXIT_FAILURE; }

    double h = 1./(double)n;    // lattice constant
    double pi   = 4.*atan(1.);  // PI
    int m  = n-1;               // matrix dimension (degrees-of-freedom)
    Matrix* diag = new Matrix(1, m, false);         // eigenvalues
    Matrix* b = new Matrix(m, m, true);             // sight side of the system
    Matrix* bt = new Matrix(m, m, true);            // tmp for transposition
    Matrix* z = new Matrix(threadMax, 4*n, false);   // tmp for FST computation

    timeStart = wallTime();

    // compute eigenvalues and right-side of the system
#pragma omp parallel for schedule(static) private(j)
    for (i=0; i < b->rows; i++) {
        diag->data[0][i] = 2.*(1.-cos((i+1)*pi/(double)n));
        for (j=0; j < b->localNumCols; j++) {
            b->data[j][i] = h*h*source((j+b->colDispl[b->worldRank]+1.0)/n,(i+1.0)/n);
        }
    }
    fstOnColumns(b, z, false);
    transpose(bt,b);;
    fstOnColumns(bt, z, true);
    // scale according to eigenvalues
#pragma omp parallel for schedule(static) private(i)
    for (j=0; j < bt->localNumCols; j++) {
        for (i=0; i < bt->rows; i++) {
            bt->data[j][i] = bt->data[j][i]/(diag->data[0][i]+diag->data[0][j+bt->colDispl[bt->worldRank]]);
        }
    }
    fstOnColumns(bt, z, false);
    transpose(b,bt);
    fstOnColumns(b, z, true);

    umax = maxError(b, &exact);
    timeEnd = wallTime();

    if (masterNode()) {
        cout << umax << "\t" << (timeEnd-timeStart) << endl;
    }
    finalize();
    return EXIT_SUCCESS;
}

void transpose (Matrix* bt, Matrix* b) {
    int row, column, elem;
    if (b->worldSize > 1) {
        double* sendbuf = new double[b->rows*b->localNumCols];
        double* recvbuf = new double[b->rows*b->localNumCols];
#pragma omp parallel for schedule(static) private(row)
        for (column = 0; column < b->localNumCols; column++) {
            for (row = 0; row < b->worldSize; row++) {
                memcpy(sendbuf+b->sendDispl[row]+column*b->numCols[row], b->data[column]+b->colDispl[row], b->numCols[row]*sizeof(double));
            }
        }
#ifdef HAVE_MPI
#pragma omp master
        MPI::COMM_WORLD.Alltoallv(sendbuf, b->sendCounts, b->sendDispl, MPI::DOUBLE, recvbuf, b->sendCounts, b->sendDispl, MPI::DOUBLE);
#endif
#pragma omp parallel for schedule(static) private(row)
        for (row = 0; row < b->worldSize; row++) {
            for (column = 0; column < b->numCols[row]; column++) {
                for (elem = 0; elem < b->localNumCols; elem++) {
                    bt->data[elem][b->colDispl[row]+column] = recvbuf[bt->sendDispl[row]+column*bt->localNumCols+elem];
                }
            }
        }
        delete[] sendbuf;
        delete[] recvbuf;
    } else {
#pragma omp for schedule(static) private(row)
        for (column=0; column < b->cols; column++) {
            for (row=0; row < b->rows; row++) {
                bt->data[column][row] = b->data[row][column];
            }
        }
    }
}

void fstOnColumns(Matrix* b, Matrix* z, bool inverse) {
    void (*fstFunction)(double *v, int *n, double *w, int *nn);
    int threadNum, j;
    int problemSize = b->rows+1;
    if (inverse) { fstFunction = &(fstinv_);
    } else { fstFunction = &(fst_); }
#pragma omp parallel private(threadNum)
    {
#ifdef HAVE_OPENMP
        threadNum = omp_get_thread_num();
#else
        threadNum = 0;
#endif
#pragma omp for schedule(static)
        for (j = 0; j < b->localNumCols; j++) {
            fstFunction(b->data[j], &problemSize, z->data[threadNum], &(z->rows));
        }
    }
}

double maxError(Matrix *b, double(*exactSolution)(double, double)) {
    double umax = 0.0;
    double tmp;
    for (int j=0; j < b->localNumCols; j++) {
        for (int i=0; i < b->rows; i++) {
            tmp = fabs(b->data[j][i]-exactSolution((j+b->colDispl[b->worldRank]+1.0)/(b->rows+1),(i+1.0)/(b->rows+1)));
            if (tmp > umax) umax = tmp;
        }
    }
#ifdef HAVE_MPI
    double tmpUmax = umax;
    MPI::COMM_WORLD.Reduce(&tmpUmax, &umax, 1, MPI::DOUBLE, MPI::MAX, 0);
#endif
    return umax;
}

bool masterNode() {
#ifdef HAVE_MPI
    return (MPI::COMM_WORLD.Get_rank() == 0);
#else
    return true;
#endif
}

bool parseArguments(int argc, char** argv, int& n) {
    if (argc < 2) {
        cerr << "error: Need a problem size." << endl;
        return false;
    }
    istringstream ss(argv[1]);
    if (!(ss >> n)) {
        cerr << "error: Could not parse problem size from '" << argv[1] << "'." << endl;
        return false;
    }
    if (n <= 0 || ( n&(n-1) != 0)) {
        cerr << "error: Problem size of " << n << " is unacceptable (need positive power of 2)." << endl;
        return false;
    }
    return true;
}

void initialize(int argc, char** argv) {
    int worldSize;
#ifdef HAVE_MPI
    MPI::Init(argc, argv);
    worldSize = MPI::COMM_WORLD.Get_size();
#else
    worldSize = 1;
#endif
#ifdef HAVE_OPENMP
    threadMax = omp_get_max_threads();
#else
    threadMax = 1;
#endif
    if (masterNode()) {
        cout << worldSize << "\t" << threadMax << "\t";
    }
}

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
