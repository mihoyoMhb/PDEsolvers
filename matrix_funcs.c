#include "matrix_funcs.h"
#define EPS 1e-6


/* 动态分配矩阵 */
double** allocate_matrix(int n, int m) {
    double **mat = (double **)malloc(n * sizeof(double *));
    if (mat == NULL) {
        perror("Matrix memory allocation failed");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n; i++) {
        mat[i] = (double *)malloc(m * sizeof(double));
        if (mat[i] == NULL) {
            perror("Matrix memory allocation failed");
            exit(EXIT_FAILURE);
        }
    }
    return mat;
}

/* freee */
void free_matrix(double **mat, int n) {
    for (int i = 0; i < n; i++) {
        free(mat[i]);
    }
    free(mat);
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

double condtion_number(double** A, int n){
    double** inv = allocate_matrix(n, n);
    inverseMatrix(A, inv, n);
    double norm_A = matrix_inf_norm(A, n);
    double norm_inv = matrix_inf_norm(inv, n);
    double cond = norm_A * norm_inv;
    free_matrix(inv, n);
    return cond;
}