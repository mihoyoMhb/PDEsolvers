#ifndef MATRIX_FUNCS_H
#define MATRIX_FUNCS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int inverseMatrix(double **A, double **inv, int n);
double matrix_inf_norm(double** A, int n);
double condtion_number(double** A, int n);
/* 动态分配 n 行 m 列的矩阵 */
double** allocate_matrix(int n, int m);

/* 释放矩阵内存 */
void free_matrix(double **mat, int n);


#endif // MATRIX_FUNCS_H
