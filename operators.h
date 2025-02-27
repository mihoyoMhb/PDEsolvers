#ifndef OPERATORS_H
#define OPERATORS_H


#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// Digaonal marix for H and HI
typedef struct{
    double *diag; // Diagonal elements
    int size; // Size of the matrix
} DiagMatrix;


// Spares matrix with CSR format (for D1 and D2)
typedef struct{
    double *values; // Non-zero values
    int *col_indices; // Column indices
    int *row_ptr; // Row pointers
    int nnz; // Number of non-zero elements
    int rows; // Number of rows
    int cols; // Number of columns
}CSRMatrix;

// Boundary vectors(e_l, e_r, d1_l, d1_r)
typedef struct{
    double* data;
    int size;
}BoundayVector;

void sbp_cent_2nd(int m, double h, 
    DiagMatrix* H, DiagMatrix* HI, 
    CSRMatrix* D1, CSRMatrix* D2,
    BoundayVector* e_l, BoundayVector* e_r,
    BoundayVector* d1_l, BoundayVector* d1_r);
#endif // OPERATORS_H