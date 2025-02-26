#include "matrix_decomposition_dynamic.h"

/* LU  */
void LUdecomposition(double **A, double **L, double **U, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            L[i][j] = (i == j) ? 1.0 : 0.0;
            U[i][j] = 0.0;
        }
    }
    for (int i = 0; i < n; i++) {
        // 计算 U 的第 i 行
        for (int j = i; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < i; k++)
                sum += L[i][k] * U[k][j];
            U[i][j] = A[i][j] - sum;
        }
        // 计算 L 的第 i 列
        for (int j = i + 1; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < i; k++)
                sum += L[j][k] * U[k][i];
            L[j][i] = (A[j][i] - sum) / U[i][i];
        }
    }
}

/* Cholesky 分解函数 */
int CholeskyDecomposition(double **A, double **L, int n) {
    // 初始化 L 为 0 矩阵
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            L[i][j] = 0.0;
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += L[i][k] * L[j][k];
            }
            if (i == j) {
                double diff = A[i][i] - sum;
                if (diff <= 0.0) {
                    return -1;  // 非正定矩阵
                }
                L[i][j] = sqrt(diff);
            } else {
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
    }
    return 0;
}

/* 打印  */
void printMatrix(double **M, int n, int m) {
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            printf("%8.4f ", M[i][j]);
        }
        printf("\n");
    }
}

/* 计算 A = L_true * L_true^T，要求 L_true 为下三角矩阵 */
void multiplyLowerTriangular(double **L_true, double **A, int n) {
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            double sum = 0.0;
            int minIndex = (i < j) ? i : j;
            for (int k = 0; k <= minIndex; k++){
                sum += L_true[i][k] * L_true[j][k];
            }
            A[i][j] = sum;
        }
    }
}
