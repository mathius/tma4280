#include <iostream>
#include <sstream>
#include <cstdlib>
#include <sys/time.h>
#include <mpi.h>

using namespace std;

#define N 5

typedef struct Vector {
    int globalLength;
    int localLength;
    int localOffset;
    int commSize;
    int commRank;
    double* data;

    Vector(int length) {
        globalLength = length;
        localLength = length;
        localOffset = 0;
        commRank = 0;
        commSize = 1;
        data = new double[localLength];
        for (int i = 0; i < localLength; i++) {
            data[i] = i;
        }

    }
    Vector(int length, int size, int rank) {
        globalLength = length;
        commSize = size;
        commRank = rank;
        localLength = globalLength / commSize;
        if (commRank < globalLength % commSize) {
            localLength++;
        }
        localOffset = commRank * (globalLength / commSize) + min(globalLength % commSize, commRank);
        data = new double[localLength];
        for (int i = 0; i < localLength; i++) {
            data[i] = i;
        }
    }
    ~Vector() {
        delete[] data;
    }
} Vector;

typedef struct Matrix {
    int nRows;
    int nColumns;
    int localNColumns;
    int localColOffset;
    double* data;
    int commSize;
    int commRank;
    Matrix(int rows, int columns, int size, int rank) {
        commSize = size;
        commRank = rank;
        nRows = rows;
        nColumns = columns;
        localNColumns = nColumns / commSize;
        if (commRank < nColumns % commSize) {
            localNColumns++;
        }
        localColOffset = commRank * (nColumns / commSize) + min(nColumns % commSize, commRank);
        data = new double[localNColumns*nRows];
        for (int i = 0; i < localNColumns*nRows; i++) {
            data[i] = i;
        }
    }
    ~Matrix() {
        delete[] data;
    }
} Matrix;

double wallTime ()
{
    struct timeval tmpTime; gettimeofday(&tmpTime,NULL);
    return tmpTime.tv_sec + tmpTime.tv_usec/1.0e6;
}

int main(int argc, char** argv) {
    srand(wallTime());
    int rank, size;

    MPI::Init(argc, argv);
    MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
    size = MPI::COMM_WORLD.Get_size();
    rank = MPI::COMM_WORLD.Get_rank();

    if (argc < 2) {
        cerr << "error: Need one numerical argument." << endl;
        return EXIT_FAILURE;
    }

    double gamma;
    stringstream ss(argv[1]);
    ss >> gamma;
    if (ss.fail()) {
        cerr << "error: Could not interpret argument '" << string(argv[1]) << "' as number." << endl;
        return EXIT_FAILURE;
    }

    Vector* a = new Vector(N, size, rank);
    Vector* b = new Vector(N, size, rank);
    Vector* x = new Vector(N, size, rank);
    Vector* y = new Vector(N);
    Vector* yTemp = new Vector(N);
    Matrix* A = new Matrix(N, N, size, rank);

    for (int i = 0; i < b->localLength; i++) {
        x->data[i] = a->data[i] + gamma * b->data[i];
        //result += x->data[i] * x->data[i];
        //cout << "rank: " << rank << ", result: " << x->data[i] << endl;
    }
    //mpiTmp = result;
    //MPI::COMM_WORLD.Allreduce(&mpiTmp, &result, 1, MPI_DOUBLE, MPI_SUM);

    for (int row = 0; row < A->nRows; row++) {
        yTemp->data[row] = 0;
        for (int i = 0; i < b->localLength; i++) {
            yTemp->data[row] += A->data[row + i * A->nRows] * b->data[i];
            //cout << "### rank: " << rank << " " << A->data[row + i * A->nRows] * b->data[i] << endl;
        }
        //cout << "rank: " << rank << ", yTemp[" << row << "]: " << yTemp->data[row] << endl;
    }
    for (int i = 0; i < yTemp->globalLength; i++) {
        y->data[i] = 0;
        MPI::COMM_WORLD.Allreduce(&(yTemp->data[i]), &(y->data[i]), 1, MPI::DOUBLE, MPI::SUM);
        //cout << "rank: " << rank << ", y[" << i << "]: " << y->data[i] << endl;
    }

    double alpha = 0;
    double alphaTmp = 0;
    for (int i = 0; i < x->localLength; i++) {
        //cout << "rank: " << rank << " alpha[" << i << "]: " << x->data[i] * (a->data[i] + y->data[i + x->localOffset]) << endl;
        alphaTmp += x->data[i] * (a->data[i] + y->data[i + x->localOffset]);
    }
    MPI::COMM_WORLD.Reduce(&alphaTmp, &alpha, 1, MPI::DOUBLE, MPI::SUM, 0);

    if (rank == 0) {
        cout << "result: " << alpha << endl;
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}

/*

a,b: 0 1 2 0 1
x = a + gamma*b = 0 1.1 2.2 0 1.1
A =
0 5 10 0 5
1 6 11 1 6
2 7 12 2 7
3 8 13 3 8
4 9 14 4 9

Ab =
(0+5+20)+(0+5) = 30
(0+6+22)+(0+6) = 34
38
42
46

y = a+Ab = 30 35 40 42 47

xTy = 0 + 1.1*35 + 2.2*40 + 0 + 1.1*47

 */