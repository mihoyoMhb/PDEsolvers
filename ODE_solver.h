#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include "matrix_funcs.h"
#include <cblas.h>
// Function pointer for ODE function

// struct with parameters for ODE function
typedef struct{
    double* v;
    double** matrix_operator;
} ODE_parameters;

// Function pointer for ODE function
typedef void (*ODE_func)(void *);
void RK45_stepping(double* t, double*v, double* dt, ODE_func f, ODE_parameters* p, int n);


#endif // ODE_SOLVER_H