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

int threadMax;

typedef struct Matrix {
    int worldSize;
    int worldRank;
    int rows;
    int cols;
    int localNumCols;
    int* numCols;
    int* colDispl;
    int* sendDispl;
    int* sendCounts;
    double* rawData;
    double** data;
    Matrix(int pCols, int pRows, bool mpiDistributed) {
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
    ~Matrix() {
        delete[] data;
        delete[] rawData;
        delete[] sendCounts;
        delete[] sendDispl;
        delete[] numCols;
    }
} Matrix;

// Sine Transform function prototypes (Fortran, use calling as C)
extern "C" {
    void fst_(double *v, int *n, double *w, int *nn);
    void fstinv_(double *v, int *n, double *w, int *nn);
}

void transpose(Matrix *bt, Matrix *b);
double maxError(Matrix* b, int n);
bool parseArguments(int argc, char** argv, int& n);
void initialize(int argc, char** argv);
void finalize();
double wallTime();

void writeMatrix(Matrix* b) {

    for (int i = 0; i < b->localNumCols; i++) {
        cout << "#" << MPI::COMM_WORLD.Get_rank() << ": ";
        for (int j = 0; j < b->rows; j++) {
            cout << b->data[i][j] << " ";
        }
        cout << endl;
    }
}

// the total number of grid points in each spatial direction is (n+1)
// the total number of degrees-of-freedom in each spatial direction is (n-1)
// this version requires n to be a power of 2

double exact(double x, double y) {
  return x*(pow(x,5)-1.0)*y*(pow(y,5)-1.0);
}

double source(double x, double y) {
  return -30.0*pow(y,4)*x*(pow(x,5.0)-1)-30.0*pow(x,4)*y*(pow(y,5)-1);
}

        double *sendbuf, *recvbuf;

int main(int argc, char **argv ) {
    // initialize MPI
    initialize(argc, argv);

    // problem size (n), matrix dimension(m), tempArray dimesion (nn), general indices (i,j)
    int i, j, n, m, nn;
    // solution (b), temp transposed solution (bt)
    // temp array for FST (z), diagonalisation vector (diag)
    Matrix *b, *bt, *z, *diag;

    double pi, h, umax, timeStart, timeEnd;
    int threadNum;

    if (!parseArguments(argc, argv, n)) {
        return EXIT_FAILURE;
    }

    m  = n-1;
    nn = 4*n;

    diag = new Matrix(1, m, false);
    b = new Matrix(m, m, true);
    bt = new Matrix(m, m, true);
    z = new Matrix(threadMax, nn, false);

    h    = 1./(double)n;
    pi   = 4.*atan(1.);

    timeStart = wallTime();

#pragma omp parallel private(threadNum)
    {
#ifdef HAVE_OPENMP
        threadNum = omp_get_thread_num();
#else
        threadNum = 0;
#endif
#pragma omp for schedule(static) private(j)
        for (i=0; i < b->rows; i++) {
            // compute eigenvalues
            diag->data[0][i] = 2.*(1.-cos((i+1)*pi/(double)n));
            // right-side of the equation (h^2*f(j,i))
            for (j=0; j < b->localNumCols; j++) {
                b->data[j][i] = h*h*source((j+b->colDispl[b->worldRank]+1.0)/n,(i+1.0)/n);
            }
        }
        // FST on columns
#pragma omp for schedule(static)
        for (j=0; j < b->localNumCols; j++) {
            fst_(b->data[j], &n, z->data[threadNum], &nn);
        }
}

        // transposition
        transpose(bt,b);
//        writeMatrix(bt);

#pragma omp parallel private(threadNum)
    {
#ifdef HAVE_OPENMP
        threadNum = omp_get_thread_num();
#else
        threadNum = 0;
#endif
        // inverse FST on rows (columns after transpose)
#pragma omp for schedule(static)
        for (i=0; i < bt->localNumCols; i++) {
            fstinv_(bt->data[i], &n, z->data[threadNum], &nn);
        }
        // scaling using eigenvalues

//writeMatrix(bt);
#pragma omp for schedule(static) private(i)
        for (j=0; j < bt->localNumCols; j++) {
            for (i=0; i < bt->rows; i++) {
                bt->data[j][i] = bt->data[j][i]/(diag->data[0][i]+diag->data[0][j+bt->colDispl[bt->worldRank]]);
            }
        }
//writeMatrix(bt);
        // FST on rows (columns after transpose)
#pragma omp for schedule(static)
        for (i=0; i < bt->localNumCols; i++) {
            fst_(bt->data[i], &n, z->data[threadNum], &nn);
        }
//writeMatrix(bt);
        }
        // transposition
        transpose(b,bt);
        // inverse FST on columns
#pragma omp parallel private(threadNum)
    {
#ifdef HAVE_OPENMP
        threadNum = omp_get_thread_num();
#else
        threadNum = 0;
#endif
//        writeMatrix(b);
#pragma omp for schedule(static)
        for (j=0; j < b->localNumCols; j++) {
            fstinv_(b->data[j], &n, z->data[threadNum], &nn);
        }
    }
//    writeMatrix(b);

    timeEnd = wallTime();
    umax = maxError(b, n);

    cout << "umax = " << umax << endl;
    cout << "time = " << (timeEnd-timeStart) << endl;

    // finalize MPI
    finalize();
    return EXIT_SUCCESS;
}

void transpose (Matrix* bt, Matrix* b) {
    int row, column, elem;
    if (b->worldSize >= 1) {
        sendbuf = new double[b->rows*b->localNumCols];
        recvbuf = new double[b->rows*b->localNumCols];
#pragma omp parallel for schedule(static) private(row)
        for (column = 0; column < b->localNumCols; column++) {
            for (row = 0; row < b->worldSize; row++) {
                memcpy(sendbuf+b->sendDispl[row]+column*b->numCols[row], b->data[column]+b->colDispl[row], b->numCols[row]*sizeof(double));
            }
        }
////        cout << "##" << MPI::COMM_WORLD.Get_rank() << ": ";
////        for (int k = 0; k < b->rows*b->localNumCols; k++) {
////            cout << sendbuf[k] << " ";
////        }
        cout << endl;
#ifdef HAVE_MPI
#pragma omp master
        MPI::COMM_WORLD.Alltoallv(sendbuf, b->sendCounts, b->sendDispl, MPI::DOUBLE, recvbuf, b->sendCounts, b->sendDispl, MPI::DOUBLE);
#endif
//        cout << "###" << MPI::COMM_WORLD.Get_rank() << ": ";
//        for (int k = 0; k < b->rows*b->localNumCols; k++) {
//            cout << recvbuf[k] << " ";
//        }
//        cout << endl;
#pragma omp parallel for schedule(static) private(row)
        for (row = 0; row < b->worldSize; row++) {
            for (column = 0; column < b->numCols[row]; column++) {
                for (elem = 0; elem < b->localNumCols; elem++) {
                    bt->data[elem][b->colDispl[row]+column] = recvbuf[bt->sendDispl[row]+column*bt->localNumCols+elem];
                }
            }
        }
    } else {
#pragma omp for schedule(static) private(row)
        for (column=0; column < b->cols; column++) {
            for (row=0; row < b->rows; row++) {
                bt->data[column][row] = b->data[row][column];
            }
        }
    }
}

double maxError(Matrix *b, int n) {
    double umax = 0.0;
    double tmp;
    for (int j=0; j < b->localNumCols; j++) {
        for (int i=0; i < b->rows; i++) {
            tmp = fabs(b->data[j][i]-exact((j+b->colDispl[b->worldRank]+1.0)/n,(i+1.0)/n));
            if (tmp > umax) umax = tmp;
        }
    }
#ifdef HAVE_MPI
    double tmpUmax = umax;
    MPI::COMM_WORLD.Reduce(&tmpUmax, &umax, 1, MPI::DOUBLE, MPI::MAX, 0);
#endif
    return umax;
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
    int myRank, worldSize;
#ifdef HAVE_MPI
    MPI::Init(argc, argv);
    myRank = MPI::COMM_WORLD.Get_rank();
    worldSize = MPI::COMM_WORLD.Get_size();
#else
    myRank = 0;
    worldSize = 1;
#endif
#ifdef HAVE_OPENMP
    threadMax = omp_get_max_threads();
#else
    threadMax = 1;
#endif
    if (myRank == 0) {
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
