#include "matrix_decomposition_dynamic.h"
#include <time.h>

int main() {
    /* LU 分解示例 */
    int nLU = 3;  // LU 分解的矩阵尺寸

    // 动态分配 LU 分解需要的矩阵
    double **A_lu = allocate_matrix(nLU, nLU);
    double **L_lu = allocate_matrix(nLU, nLU);
    double **U_lu = allocate_matrix(nLU, nLU);

    // 示例数据
    double A_lu_data[3][3] = {
        {4, -2, 1},
        {-2, 4, -2},
        {1, -2, 3}
    };
    for (int i = 0; i < nLU; i++){
        for (int j = 0; j < nLU; j++){
            A_lu[i][j] = A_lu_data[i][j];
        }
    }

    // 调用 LU 分解函数
    LUdecomposition(A_lu, L_lu, U_lu, nLU);

    printf("LU 分解结果：\n");
    printf("L 矩阵：\n");
    printMatrix(L_lu, nLU, nLU);
    printf("\nU 矩阵：\n");
    printMatrix(U_lu, nLU, nLU);

    /* Cholesky 分解示例 */
    int nChol = 5;  // Cholesky 分解的矩阵尺寸

    // 动态分配 Cholesky 分解需要的矩阵
    double **L_true = allocate_matrix(nChol, nChol);
    double **A_chol = allocate_matrix(nChol, nChol);
    double **L_chol = allocate_matrix(nChol, nChol);

    // 构造下三角矩阵 L_true（随机整数填充，范围 0~9）
    srand((unsigned)time(NULL));
    for (int i = 0; i < nChol; i++){
        for (int j = 0; j < nChol; j++){
            if (j <= i)
                L_true[i][j] = rand() % 10;
            else
                L_true[i][j] = 0.0;
        }
    }

    // 计算 A_chol = L_true * L_true^T
    multiplyLowerTriangular(L_true, A_chol, nChol);

    // 调用 Cholesky 分解函数
    if (CholeskyDecomposition(A_chol, L_chol, nChol) != 0){
        printf("\nCholesky 分解失败：矩阵非正定！\n");
    } else {
        printf("\nCholesky 分解结果：\n");
        printf("原始下三角矩阵 L_true：\n");
        printMatrix(L_true, nChol, nChol);
        printf("\n由 L_true * L_true^T 得到的矩阵 A：\n");
        printMatrix(A_chol, nChol, nChol);
        printf("\n从 A 分解得到的下三角矩阵 L：\n");
        printMatrix(L_chol, nChol, nChol);
    }

    /* 释放所有内存 */
    free_matrix(A_lu, nLU);
    free_matrix(L_lu, nLU);
    free_matrix(U_lu, nLU);
    free_matrix(L_true, nChol);
    free_matrix(A_chol, nChol);
    free_matrix(L_chol, nChol);

    return 0;
}
