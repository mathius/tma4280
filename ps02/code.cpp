#include <iostream>
#include <sstream>
#include <cstdlib>
#include <sys/time.h>

using namespace std;

#define N 500

double wallTime ()
{
    struct timeval tmpTime; gettimeofday(&tmpTime,NULL);
    return tmpTime.tv_sec + tmpTime.tv_usec/1.0e6;
}

int main(int argc, char** argv) {
    double* x = (double*) malloc(N*sizeof(double));
    double* y = (double*) malloc(N*sizeof(double));
    double* a = (double*) malloc(N*sizeof(double));
    double* b = (double*) malloc(N*sizeof(double));
    double* matrix = (double*) malloc(N*N*sizeof(double));
    double gamma = rand()/RAND_MAX;
    double alpha;
    int i,j;
    double before, after;

    if (x == NULL || y == NULL || a == NULL || b == NULL || matrix == NULL) {
        cerr << "error: Memory allocation failed." << endl;
        return EXIT_FAILURE;
    }

    if (argc < 2) {
        cerr << "error: Need one numerical argument." << endl;
        return EXIT_FAILURE;
    }
    stringstream ss(argv[1]);
    ss >> gamma;
    if (ss.fail()) {
        cerr << "error: Could not interpret argument '" << string(argv[1]) << "' as number." << endl;
        return EXIT_FAILURE;
    }

    srand(wallTime());

    // initialize numbers
    for (i = 0; i < N; i++) {
        a[i] = (double) rand()/RAND_MAX;
        b[i] = (double) rand()/RAND_MAX;
        for (j = 0; j < N; j++) {
            matrix[i*N+j] = (double) rand()/RAND_MAX;
        }
    }

    before = wallTime();

    // compute gamma*b to x, matrix*b to y
    for (i = 0; i < N; i++) {
        x[i] = gamma*b[i];
    }
    for (i = 0; i < N; i++) {
        y[i] = 0;
        for (j = 0; j < N; j++) {
            y[i] += matrix[N*i+j]*b[j];
        }
    }

    // add a to both x and y
    for (i = 0; i < N; i++) {
        x[i] += a[i];
        y[i] += b[i];
    }

    // compute aplha
    alpha = 0;
    for (i = 0; i < N; i++) {
        alpha += x[i]*y[i];
    }

    after = wallTime();

    cout << "result: " << alpha << endl;
    cout << "duration: " << (after-before) << endl;

    free(x);
    free(y);
    free(matrix);
    free(a);
    free(b);

    return EXIT_SUCCESS;
}