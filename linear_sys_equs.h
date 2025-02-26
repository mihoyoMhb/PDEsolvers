#ifndef LINEAR_SYS_EQUS_H
#define LINEAR_SYS_EQUS_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*Solve Ly = b*/
void lower_TriangularSolver(double **L, double *b, double* y, int n);
/*Solve L^T x = y */
void Cholesky_upper_triangularSolver(double **L, double *y, double*x, int n);
/*Solve Ax = b */
void Cholesky_Solver(double **A, double* b, double* x, int n);
/*Solve Ux=y*/
void normal_UpperTriangularSolver(double** U, double* b, double* x, int n);
#endif // LINEAR_SYS_EQUS_H