#ifndef _POISSON_H_
#define _POISSON_H_

/*
 * TMA4280 Supercomputing, Introduction
 * Problem set 6
 *
 * for authors and more info see poisson.cpp
 */

typedef struct Matrix {
    int worldSize;      //! number of nodes the matrix is distributed across
    int worldRank;      //! rank of this node
    int rows;           //! global number of row
    int cols;           //! global number of columns
    int localNumCols;   //! local number of columns
    int* numCols;       //! number of columns across nodes
    int* colDispl;      //! column displacements across nodes
    int* sendDispl;     //! send displacements for the transpose operation
    int* sendCounts;    //! send counts for the transpose operation
    double* rawData;    //! actual data, column-wise memory layout
    double** data;      //! individual columns access
    /**
     * @brief Matrix constructor
     * @param pCols     global number of columns
     * @param pRows     global number of rows
     * @param mpiDistributed    should this matrix be MPI distributed (if supported)?
     */
    Matrix(int pCols, int pRows, bool mpiDistributed);
    /**
      * @brief Matrix destructor, frees memory
      */
    ~Matrix();
} Matrix;

/**
 * @brief source function f, must have homogenous Dirichlet boundary conditions
 *        i.e. f(x,y)=0 for (x=0 || y=0)
 * @param x, y
 */
double source(double x, double y) {
  return -30.0*pow(y,4)*x*(pow(x,5.0)-1)-30.0*pow(x,4)*y*(pow(y,5)-1);
}

/**
 * @brief exact solution (for computing absolute error)
 * @param x, y
 */
double exact(double x, double y) {
  return x*(pow(x,5)-1.0)*y*(pow(y,5)-1.0);
}

// sine transform function prototypes (Fortran, use calling as C)
extern "C" {
    /**
     * @brief Fast Sine Transform
     * @param v     vector for transformation
     * @param n     problem size (vector size + 1)
     * @param w     temporary array
     * @param nn    length of temprary array (4*problemSize)
     */
    void fst_(double *v, int *n, double *w, int *nn);

    /**
     * @brief Inverse Fast Sine Transform
     * @param v     vector for transformation
     * @param n     problem size (vector size + 1)
     * @param w     temporary array
     * @param nn    length of temprary array (4*problemSize)
     */
    void fstinv_(double *v, int *n, double *w, int *nn);
}

/**
 * @brief (Inverse) Fast Sine Transform on matrix columns
 * @param b         input matrix of (problemSize - 1) rows
 * @param z         matrix of temprary arrays (4*problemSize x numberOfThreads)
 * @param inverse   should we do inverse transform?
 */
void fstOnColumns(Matrix* b, Matrix* z, bool inverse);

/**
 * @brief matrix transposition
 * @param bt    transposed matrix
 * @param b     matrix to transpose
 */
void transpose(Matrix *bt, Matrix *b);

/**
 * @brief determine maximum absolute error from exact solution
 * @param b                 matrix with computed solution
 * @param exactSolution     callback for exact solution
 * @return max absolute error
 */
double maxError(Matrix* b, double(*exactSolution)(double, double));

/**
 * @brief check the second argument is a positive power of 2 and parse it to n
 * @param argc  num of arguments
 * @param argv  cli arguments
 * @param n     parsed problem size
 * @return      parsing process success
 */
bool parseArguments(int argc, char** argv, int& n);

/**
 * @brief initialize MPI/OpenMP, fill global variables, print debug info
 * @param argc  num of arguments
 * @param argv  cli arguments
 */
void initialize(int argc, char** argv);

/**
 * @brief finalize MPI
 */
void finalize();

/**
 * @brief inspects MPI world rank
 * @return is the current node master?
 */
bool masterNode();

/**
 * @brief timing utility (MPI/OpenMP compatible)
 * @return wall time in seconds
 */
double wallTime();

#endif
