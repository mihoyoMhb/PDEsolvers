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
}BoundaryVector;

void sbp_cent_2nd(int m, double h, 
    DiagMatrix* H, DiagMatrix* HI, 
    CSRMatrix* D1, CSRMatrix* D2,
    BoundaryVector* e_l, BoundaryVector* e_r,
    BoundaryVector* d1_l, BoundaryVector* d1_r);


// CSR matrix-vector product
void csr_matvec(CSRMatrix* A, double* x, double* y);

// Diagonal matrix-vector product
void diag_matvec(DiagMatrix* A, double* x, double* y);

// Diagonal matrix-Boundary vector product
// return a boundary vector y = A*b
void diag_bv_product(DiagMatrix* A, BoundaryVector* b, BoundaryVector* y);

// Boundary vector-Boundary vector product
// For FDM, return a CSR matrix A = B1*B2
void bv_bv_product(BoundaryVector* B1, BoundaryVector* B2, CSRMatrix* A);

// CSR matrix add operation
// return a CSR matrix C = A + B
void csr_add(CSRMatrix* A, CSRMatrix* B, CSRMatrix* C);


// Free memory
void free_diag(DiagMatrix *mat);

void free_csr(CSRMatrix *mat); 

void free_boundary(BoundaryVector *vec); 
#endif // OPERATORS_H