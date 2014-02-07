#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <cstdlib>
#include <cmath>

typedef struct Vector {
    int lengthGlobal;
    int lengthLocal;
    int offsetLocal;
    int commSize;
    int commRank;
    double* data;
    Vector(int length, double* existingData = NULL);
    ~Vector();
} Vector;

void initialize(int argc, char** argv);
void finalize();

int getWorldRank();
int getWorldSize();

//! \brief Get current wall-clock time
//! \return The current wall time in seconds
double wallTime();

#endif // ALGEBRA_H
