#include "operators.h"

void csr_matvec(CSRMatrix* A, double* x, double* y){
    /*
    Compute matrix-vector product y = A*x
    A: Sparse matrix in CSR format
    x: Input vector
    y: Output vector
    */
    for(int i=0; i<A->rows; i++){
        y[i] = 0.0;
        for(int j = A->row_ptr[i]; j< A->row_ptr[i+1]; j++){
            y[i] += A->values[j] * x[A->col_indices[j]];
        }
    }
}

void diag_matvec(DiagMatrix* A, double *x, double *y){
    /*
    Compute matrix-vector product y = A*x
    A: Diagonal matrix
    x: Input vector
    y: Output vector
    */
    for(int i=0;i<A->size;i++){
        y[i] = A->diag[i] * x[i];
    }
}

// Free memory
void free_diag(DiagMatrix *mat) {
    free(mat->diag);
    mat->size = 0;
}

void free_csr(CSRMatrix *mat) {
    free(mat->values);
    free(mat->col_indices);
    free(mat->row_ptr);
    mat->rows = mat->cols = mat->nnz = 0;
}

void free_boundary(BoundaryVector *vec) {
    free(vec->data);
    vec->size = 0;
}


void sbp_cent_2nd(int m, double h, 
    DiagMatrix* H, DiagMatrix* HI, 
    CSRMatrix* D1, CSRMatrix* D2,
    BoundaryVector* e_l, BoundaryVector* e_r,
    BoundaryVector* d1_l, BoundaryVector* d1_r) {
    /*
    Construct SBP operators for 2nd order central finite difference
    m: number of grid points
    h: grid spacing
    H: Diagonal matrix for H
    HI: Diagonal matrix for HI
    D1: Sparse matrix for D1
    D2: Sparse matrix for D2
    e_l: Boundary vector e_l
    e_r: Boundary vector e_r
    d1_l: Boundary vector d1_l
    d1_r: Boundary vector d1_r
    */
    // Allocate memory for H and HI
    H->diag = (double*)malloc(m * sizeof(double));
    H->size = m;
    HI->diag = (double*)malloc(m * sizeof(double));
    HI->size = m;
    // Initialize H and HI
    for(int i=0;i<m;i++) H->diag[i] = (i == 0 || i == m-1) ? 0.5*h : h;
    for(int i=0;i<m;i++) HI->diag[i] = 1.0 / H->diag[i];
    
    // Initialize D1 with CSR format
    int nnz = 2*m;
    D1->values = (double*)malloc(nnz * sizeof(double));
    D1->col_indices = (int*)malloc(nnz * sizeof(int));
    D1->row_ptr = (int*)malloc((m+1) * sizeof(int));

    int idx = 0;
    D1->row_ptr[0] = 0;
    for(int i=0;i<m;i++){
        if(i == 0){
            D1->values[idx] = -1.0 / h;
            D1->col_indices[idx] = 0;
            idx++;
            D1->values[idx] = 1.0 / h;
            D1->col_indices[idx] = 1;
            idx++;
        }else if (i == m-1){
            D1->values[idx] = -1.0 / h;
            D1->col_indices[idx] = m-2;
            idx++;
            D1->values[idx] = 1.0 / h;
            D1->col_indices[idx] = m-1;
            idx++;
        }else{
            D1->values[idx] = -0.5 / h;
            D1->col_indices[idx] = i-1;
            idx++;
            D1->values[idx] = 0.5 / h;
            D1->col_indices[idx] = i+1;
            idx++;
        }
        D1->row_ptr[i+1] = idx;
    }
    D1->nnz = nnz;
    D1->rows = m;
    D1->cols = m;

    // Initialize e_l, e_r
    e_l->data = (double*)malloc(m * sizeof(double)); e_l->data[0] = 1.0; e_l->size = m;
    e_r->data = (double*)malloc(m * sizeof(double)); e_r->data[m-1] = 1.0; e_r->size = m;

    // Initialize d1_l, d1_r
    d1_l->data = (double*)malloc(m * sizeof(double)); d1_l->size = m;
    d1_r->data = (double*)malloc(m * sizeof(double)); d1_r->size = m;
    d1_l->data[0] = -1.5; d1_r->data[m-1] = 1.5;
    d1_l->data[1] = 2.0; d1_r->data[m-2] = -2.0;
    d1_l->data[2] = -0.5; d1_r->data[m-3] = 0.5;
    for(int i=3;i<m;i++){
        d1_l->data[i] = 0.0;
        d1_r->data[m-i-1] = 0.0;
    }

    // Initialize D2 with CSR format
    nnz = 3*m;
    D2->values = (double*)malloc(nnz * sizeof(double));
    D2->col_indices = (int*)malloc(nnz * sizeof(int));
    D2->row_ptr = (int*)malloc((m+1) * sizeof(int));
    idx = 0;
    double h_2 = h * h;
    D2->row_ptr[0] = 0;
    for(int i = 0; i < m; i++){
        if(i == 0){
            D2->values[idx] = 1.0 / h_2;
            D2->col_indices[idx] = 0;
            idx++;
            D2->values[idx] = -2.0 / h_2;
            D2->col_indices[idx] = 1;
            idx++;
            D2->values[idx] = 1.0 / h_2;
            D2->col_indices[idx] = 2;
            idx++;
        }else if(i == m-1){
            D2->values[idx] = 1.0 / h_2;
            D2->col_indices[idx] = m-3;
            idx++;
            D2->values[idx] = -2.0 / h_2;
            D2->col_indices[idx] = m-2;
            idx++;
            D2->values[idx] = 1.0 / h_2;
            D2->col_indices[idx] = m-1;
            idx++;
        }else{
            D2->values[idx] = 1.0 / h_2;
            D2->col_indices[idx] = i-1;
            idx++;
            D2->values[idx] = -2.0 / h_2;
            D2->col_indices[idx] = i;
            idx++;
            D2->values[idx] = 1.0 / h_2;
            D2->col_indices[idx] = i+1;
            idx++;
        }
        D2->row_ptr[i+1] = idx;
    }
    D2->nnz = nnz;
    D2->rows = m;
    D2->cols = m;
}