#include "operators.h"

void sbp_cent_2nd(int m, double h, 
    DiagMatrix* H, DiagMatrix* HI, 
    CSRMatrix* D1, CSRMatrix* D2,
    BoundayVector* e_l, BoundayVector* e_r,
    BoundayVector* d1_l, BoundayVector* d1_r) {
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
    

}