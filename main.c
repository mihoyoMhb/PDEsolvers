#include "matrix_decomposition_dynamic.h"
#include "linear_sys_equs.h"
#include "matrix_funcs.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h> 
#define EPS 1e-6

/* 辅助函数：计算矩阵与向量的乘积，r = A*x */
void matVecProduct(double **A, double *x, double *r, int n) {
    for (int i = 0; i < n; i++) {
        r[i] = 0.0;
        for (int j = 0; j < n; j++) {
            r[i] += A[i][j] * x[j];
        }
    }
}


void printVector(double *v, int n) {
    for (int i = 0; i < n; i++) {
        printf("%8.4f ", v[i]);
    }
    printf("\n");
}


int main(int argc, char const *argv[]) {
    if(argc != 2){
        printf("Usage: %s <size of matrix>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);  // 
    
    srand((unsigned)time(NULL));

    double **A = allocate_matrix(n, n);
    double **L_true = allocate_matrix(n, n);

    /* 构造随机下三角矩阵 L_true */
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (j <= i) {
                // 生成 -1 到 1 之间的随机 double 数值
                double rnd = 2.0 * rand() / RAND_MAX - 1.0;
                // 如果是对角元，可以确保其不为 0（例如重新生成直到不为 0）
                if (i == j) {
                    while(rnd == 0.0) {
                        rnd = 2.0 * rand() / RAND_MAX - 1.0;
                    }
                }
                L_true[i][j] = rnd;
            } else {
                L_true[i][j] = 0.0;
            }
        }
    }

    /* 计算 A = L_true * L_true^T */
    multiplyLowerTriangular(L_true, A, n);

    // printf("生成的对称正定矩阵 A:\n");
    // printMatrix(A, n, n);
    // printf("\n");

    /* 随机生成右端向量 b */
    double *b = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++){
        b[i] = 2.0 * rand() / RAND_MAX - 1.0;  // 取 -1~1 内整数
    }
    // printf("右端向量 b:\n");
    // printVector(b, n);
    // printf("\n");

    printf("Condition number of A: %f\n", condtion_number(A, n));
    /* 调用 Cholesky_Solver 求解 Ax = b */
    double *x = (double*)malloc(n * sizeof(double));
    Cholesky_Solver(A, b, x, n);

    printf("求解得到的解向量 x:\n");
    printVector(x, n);
    printf("\n");

    /* 验证解：计算 r = A*x - b */
    double *r = (double*)malloc(n * sizeof(double));
    matVecProduct(A, x, r, n);
    for (int i = 0; i < n; i++){
        r[i] -= b[i];
    }
    printf("残差 (A*x - b):\n");
    printVector(r, n);
    printf("\n");

    /* 释放内存 */
    free_matrix(A, n);
    free_matrix(L_true, n);
    free(b);
    free(x);
    free(r);

    return 0;
}
