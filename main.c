#include <stdio.h>
#include <stdlib.h>

#define N 3  // 假设矩阵是 N x N

void LUdecomposition(double A[N][N], double L[N][N], double U[N][N]) {
    for (int i = 0; i < N; i++) {
        // 计算 U 的第 i 行
        for (int j = i; j < N; j++) {
            double sum = 0;
            for (int k = 0; k < i; k++)
                sum += L[i][k] * U[k][j];
            U[i][j] = A[i][j] - sum;
        }

        // 计算 L 的第 i 列
        for (int j = i + 1; j < N; j++) {
            double sum = 0;
            for (int k = 0; k < i; k++)
                sum += L[j][k] * U[k][i];
            L[j][i] = (A[j][i] - sum) / U[i][i];
        }
    }
}

void printMatrix(double M[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            printf("%8.4f ", M[i][j]);
        printf("\n");
    }
}

int main() {
    double A[N][N] = {
        {4, -2, 1},
        {-2, 4, -2},
        {1, -2, 3}
    };

    double L[N][N] = {0}, U[N][N] = {0};

    // 初始化 L 为单位矩阵
    for (int i = 0; i < N; i++)
        L[i][i] = 1;

    LUdecomposition(A, L, U);

    printf("L matrix:\n");
    printMatrix(L);
    printf("\nU matrix:\n");
    printMatrix(U);

    return 0;
}
