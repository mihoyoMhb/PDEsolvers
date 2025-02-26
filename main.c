#include "matrix_decomposition_dynamic.h"
#include "linear_sys_equs.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* 辅助函数：计算矩阵与向量的乘积，r = A*x */
void matVecProduct(double **A, double *x, double *r, int n) {
    for (int i = 0; i < n; i++) {
        r[i] = 0.0;
        for (int j = 0; j < n; j++) {
            r[i] += A[i][j] * x[j];
        }
    }
}

/* 辅助函数：打印向量 */
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
    int n = atoi(argv[1]);  // 矩阵尺寸
    
    srand((unsigned)time(NULL));

    /* 动态分配矩阵 A 和下三角矩阵 L_true */
    double **A = allocate_matrix(n, n);
    double **L_true = allocate_matrix(n, n);

    /* 构造随机下三角矩阵 L_true
       为保证正定性，设定对角元在 [1,10] 内，其余下三角元在 [0,10] 内 */
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (j <= i) {
                if(i == j)
                    L_true[i][j] = 1 + rand() % 10;  // 对角元不为 0
                else
                    L_true[i][j] = rand() % 10;
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
        b[i] = rand() % 10;  // 取 0~9 内整数
    }
    // printf("右端向量 b:\n");
    // printVector(b, n);
    // printf("\n");

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
