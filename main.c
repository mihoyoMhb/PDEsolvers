#include "matrix_decomposition_dynamic.h"
#include "linear_sys_equs.h"
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


/* 计算矩阵 A 的逆矩阵，结果存入 inv 中，A 与 inv 均为 n×n 的动态矩阵
   若矩阵 A 为奇异矩阵，则返回 -1 */
int inverseMatrix(double **A, double **inv, int n) {
    int i, j, k;
    double pivot, ratio;
    // 分配增广矩阵 aug，尺寸为 n x (2*n)
    double **aug = allocate_matrix(n, 2 * n);
    if (!aug) {
        printf("无法分配增广矩阵！\n");
        return -1;
    }
    // 构造增广矩阵 [A | I]
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            aug[i][j] = A[i][j];
        }
        for (j = n; j < 2 * n; j++) {
            aug[i][j] = (i == (j - n)) ? 1.0 : 0.0;
        }
    }
    // 高斯－约旦消元
    for (i = 0; i < n; i++) {
        // 检查主元是否足够大
        if (fabs(aug[i][i]) < EPS) {
            int swapRow = -1;
            for (j = i + 1; j < n; j++) {
                if (fabs(aug[j][i]) > EPS) {
                    swapRow = j;
                    break;
                }
            }
            if (swapRow == -1) {
                free_matrix(aug, n);
                return -1;  // 矩阵奇异，无法求逆
            }
            // 交换第 i 行与 swapRow 行
            for (k = 0; k < 2 * n; k++) {
                double temp = aug[i][k];
                aug[i][k] = aug[swapRow][k];
                aug[swapRow][k] = temp;
            }
        }
        // 将第 i 行归一化，使主元为 1
        pivot = aug[i][i];
        for (j = 0; j < 2 * n; j++) {
            aug[i][j] /= pivot;
        }
        // 对其他各行消元
        for (j = 0; j < n; j++) {
            if (j != i) {
                ratio = aug[j][i];
                for (k = 0; k < 2 * n; k++) {
                    aug[j][k] -= ratio * aug[i][k];
                }
            }
        }
    }
    // 从增广矩阵中提取逆矩阵，即右半部分
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            inv[i][j] = aug[i][j + n];
        }
    }
    free_matrix(aug, n);
    return 0;
}

/*INF norm*/
double matrix_inf_norm(double** A, int n){
    double row_sum[n];
    for(int i=0; i<n; i++){
        row_sum[i] = 0;
        for(int j=0; j<n; j++){
            row_sum[i] += abs(A[i][j]);
        }
    }
    double max_row_sum = -INFINITY;
    for(int i=1; i<n; i++){
        if(row_sum[i] > max_row_sum){
            max_row_sum = row_sum[i];
        }
    }
    return max_row_sum;
}

void condtion_number(double** A, int n){
    double** inv = allocate_matrix(n, n);
    inverseMatrix(A, inv, n);
    double norm_A = matrix_inf_norm(A, n);
    double norm_inv = matrix_inf_norm(inv, n);
    double cond = norm_A * norm_inv;
    printf("The condition number of A is %f\n", cond);
    free_matrix(inv, n);
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

    condtion_number(A, n);

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
