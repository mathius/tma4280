#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

#define N 300

int main(int argc, char** argv) {
	srand(time(NULL));

	double* x = (double*) malloc(N*sizeof(double));
	double* y = (double*) malloc(N*sizeof(double));
	double* a = (double*) malloc(N*N*sizeof(double));

	for (int i = 0; i < N; i++) {
		x[i] = (double) rand()/RAND_MAX;
		for (int j = 0; j < N; j++) {
			a[i*N+j] = (double) rand()/RAND_MAX;
		}
	}

	unsigned long before =  time(NULL);

	for (int i = 0; i < N; i++) {
		y[i] = 0;
		for (int j = 0; j < N; j++) {
			y[i] += a[i*N+j]*x[i];
		}
	}

	unsigned long after = time(NULL);

    cout << "Duration: " << (after-before) << endl;
    cout << "Result: ";
	for (int i = 0; i < N; i++) {
		cout << y[i] << " ";
	}
	cout << endl;

	free(x);
    free(y);
    free(a);

	return EXIT_SUCCESS;
}
