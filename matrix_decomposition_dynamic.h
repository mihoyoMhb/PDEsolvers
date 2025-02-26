#ifndef MATRIX_DECOMPOSITION_DYNAMIC_H
#define MATRIX_DECOMPOSITION_DYNAMIC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* 动态分配 n 行 m 列的矩阵 */
double** allocate_matrix(int n, int m);

/* 释放矩阵内存 */
void free_matrix(double **mat, int n);


void LUdecomposition(double **A, double **L, double **U, int n);

/*
 * Cholesky 分解函数
 * 输入：A 为 n×n 正定对称矩阵，L 为输出的下三角矩阵（n×n），
 * 返回值：若矩阵非正定，返回 -1；否则返回 0。
 */
int CholeskyDecomposition(double **A, double **L, int n);

/* 打印  */
void printMatrix(double **M, int n, int m);

/*
 * 计算 A = L_true * L_true^T，
 * 其中 L_true 为 n×n 下三角矩阵。
 */
void multiplyLowerTriangular(double **L_true, double **A, int n);

#endif // MATRIX_DECOMPOSITION_DYNAMIC_H
