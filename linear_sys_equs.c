#include "linear_sys_equs.h"
#include "matrix_decomposition_dynamic.h"

/*Solve Ly=b*/
void lower_TriangularSolver(double** L, double* b, double* y, int n){
    for(int i=0; i<n; i++){
        y[i] = b[i];
        for(int j=0; j<i; j++){
            y[i] -= L[i][j]*y[j];
        }
        y[i] /= L[i][i];
    }
}

/*Solve L^T x = y */
void Cholesky_upper_triangularSolver(double** L, double* y, double* x, int n){
    /*Since A = LL^T, we can easily to calculate the function*/
    for(int i=n-1; i>=0; i--){
        x[i] = y[i];
        for(int j=i+1; j<n; j++){
            x[i] -= L[j][i]*x[j];
        }
        x[i] /= L[i][i];
    }
}

/*Solve Ux=y*/
void normal_UpperTriangularSolver(double** U, double* b, double* x, int n){
    for(int i=n-1; i>=0; i--){
        x[i] = b[i];
        for(int j=i+1; j<n; j++){
            x[i] -= U[i][j]*x[j];
        }
        x[i] /= U[i][i];
    }
}


/*Solve Ax=b*/
void Cholesky_Solver(double** A, double* b, double* x, int n){
    double y[n];
    double **L = allocate_matrix(n, n);
    CholeskyDecomposition(A, L, n);
    lower_TriangularSolver(L, b, y, n);
    Cholesky_upper_triangularSolver(L, y, x, n);
    free_matrix(L, n);
}