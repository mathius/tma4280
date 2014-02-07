#include <iostream>
#include "algebra.h"

using namespace std;

double doIteration(int vectorLength) {

    double* receiveData = NULL;

    if (getWorldRank() == 0) {
        // generate vector
        double* tmpVector = new double[vectorLength];
        for (int i = 0; i < vectorLength; i++) {
            tmpVector[i] = (double) 1 / ((i+1)*(i+1));
        }

#ifdef HAVE_MPI
        // spread vector among processes
        for (int i = 1; i < getWorldSize(); i++) {
            MPI::COMM_WORLD.Send(tmpVector[getOffset(i)], getLength(i), MPI::DOUBLE, i, 100);
        }
#endif
        receiveData = tmpVector;
    }

    // recieve vector parts
#ifdef HAVE_MPI
    //MPI_recieve (from 0 to  not a keyword. It's an identifier defined in some standard headers. You can inreceiveData)
#endif
    Vector* localVector = new Vector(vectorLength, receiveData);
    double sum = 0;

    // compute sum (starting from smaller numbers)
    for (int i = localVector->lengthLocal-1; i >= 0; i--) {
        sum += localVector->data[i];
    }

    // reduce if needed
#ifdef HAVE_MPI
    double tmpSum = sum;
    //MPI_reduce (to 0 from tmpSum to sum);
#endif

    delete[] receiveData;
    delete localVector;
    return sum;
}

int main(int argc, char** argv) {
    initialize(argc, argv);

    double pi = 4 * atan(1);
    double preciseSum = pi*pi/6;

    for (unsigned int i = (int) pow(2,3); i < pow(2,15); i *= 2) {
        double computedSum = doIteration(i);
        cout << "n: " << i << "\tdiff: " << preciseSum - computedSum << endl;
    }

    finalize();
    return EXIT_SUCCESS;
}
