#include <math.h>

#include "log.h"

#include "linalg.h"

#define MAX_ITER 100
#define EPSILON 1e-12

static void product_AtA(const matrix_t *A, matrix_t *B)
{
    size_t rows;
    size_t cols;

    rows = A->size1;
    cols = A->size2;

    for (size_t i = 0; i < cols; i++) {
        for (size_t j = i; j < cols; j++) {
            double sum;

            sum = 0.0;
            for (size_t k = 0; k < rows; k++)
                sum += A->data[k * cols + i] * A->data[k * cols + j];

            matrix_set(B, i, j, sum);
            matrix_set(B, j, i, sum);
        }
    }
}

static void jacobi_svd(matrix_t *B, matrix_t *V, vector_t *S)
{
    size_t cols;
    double max_offdiag;
    size_t iter;

    cols = B->size1;

    matrix_set_identity(V);

    for (iter = 0; iter < MAX_ITER; iter++) {
        max_offdiag = 0.0;

        for (size_t p = 0; p < cols - 1; p++) {
            for (size_t q = p + 1; q < cols; q++) {
                double B_pp;
                double B_qq;
                double B_pq;
                double tau;
                double t;
                double c;
                double s;

                B_pp = B->data[p * cols + p];
                B_qq = B->data[q * cols + q];
                B_pq = B->data[p * cols + q];

                if (fabs(B_pq) < EPSILON)
                    continue;

                tau = (B_qq - B_pp) / (2.0 * B_pq);
                if (tau >= 0.0)
                    t = 1.0 / (tau + sqrt(1.0 + tau * tau));
                else
                    t = 1.0 / (tau - sqrt(1.0 + tau * tau));

                c = 1.0 / sqrt(1.0 + t * t);
                s = t * c;

                // Update rows p and q of B
                for (size_t k = 0; k < cols; k++) {
                    double B_pk;
                    double B_qk;

                    B_pk = B->data[p * cols + k];
                    B_qk = B->data[q * cols + k];

                    B->data[p * cols + k] = c * B_pk - s * B_qk;
                    B->data[q * cols + k] = s * B_pk + c * B_qk;
                }

                // Update columns p and q of B
                for (size_t k = 0; k < cols; k++) {
                    double B_kp;
                    double B_kq;

                    B_kp = B->data[k * cols + p];
                    B_kq = B->data[k * cols + q];

                    B->data[k * cols + p] = c * B_kp - s * B_kq;
                    B->data[k * cols + q] = s * B_kp + c * B_kq;
                }

                // Update columns p and q of V
                for (size_t k = 0; k < cols; k++) {
                    double V_kp;
                    double V_kq;

                    V_kp = V->data[k * cols + p];
                    V_kq = V->data[k * cols + q];

                    V->data[k * cols + p] = c * V_kp - s * V_kq;
                    V->data[k * cols + q] = s * V_kp + c * V_kq;
                }

                if (fabs(B_pq) > max_offdiag)
                    max_offdiag = fabs(B_pq);
            }
        }

        if (max_offdiag < EPSILON)
            break;
    }

    // Extract singular values as square roots of the diagonal of B
    for (size_t i = 0; i < cols; i++) {
        double val;

        val = B->data[i * cols + i];
        if (val > 0.0)
            S->data[i] = sqrt(val);
        else
            S->data[i] = 0.0;
    }
}

static void compute_U_from_AV(const matrix_t *A_orig, const matrix_t *V, const vector_t *S, matrix_t *U_out)
{
    size_t rows;
    size_t cols;

    rows = A_orig->size1;
    cols = A_orig->size2;

    for (size_t i = 0; i < cols; i++) {
        double sigma;

        sigma = S->data[i];
        if (sigma < EPSILON) {
            // Set column i of U to zero if sigma is too small
            for (size_t row = 0; row < rows; row++)
                U_out->data[row * cols + i] = 0.0;
            continue;
        }

        for (size_t row = 0; row < rows; row++) {
            double sum;

            sum = 0.0;
            for (size_t j = 0; j < cols; j++)
                sum += A_orig->data[row * cols + j] * V->data[j * cols + i];

            U_out->data[row * cols + i] = sum / sigma;
        }
    }
}

/**
 * @brief Performs LU decomposition with partial pivoting.
 *
 * This function decomposes a matrix A into a lower trapezoidal matrix L and an
 * upper trapezoidal matrix U such that PA = LU, where P is a permutation matrix.
 * The permutation information is stored in the permutation_t structure.
 * The LU matrix is stored in place of A.
 *
 * @param A A matrix to be decomposed. On output, it contains L and U.
 * @param p The permutation object to store row interchanges.
 * @param signum A pointer to an integer that will hold the sign of the permutation.
 * @return 0 on success, -1 on failure (e.g., matrix is singular).
 */
int linalg_lu_decomp(matrix_t *A, permutation_t *p, int *signum)
{
    size_t dim_min;
    size_t rows;
    size_t cols;

    rows = A->size1;
    cols = A->size2;
    dim_min = (rows < cols) ? rows : cols;

    if (p->size < rows) {
        LOG_ERROR("Error: Permutation vector size is too small");
        return -1;
    }

    * signum = 1;
    permutation_set_identity(p);

    for (size_t k = 0; k < dim_min; k++) {
        /* Pivoting: biggest elem in column k */
        size_t row_max;
        double val_max;
        double val;

        row_max = k;
        val_max = fabs(A->data[k * cols + k]);
        for (size_t i = k + 1; i < rows; i++) {
            val = fabs(A->data[i * cols + k]);
            if (val > val_max) {
                val_max = val;
                row_max = i;
            }
        }

        /* Swap rows if pivot is not in current row */
        if (row_max != k) {
            matrix_swap_rows(A, k, row_max);
            permutation_swap(p, k, row_max);
            *signum = -(*signum);
        }

        /* Gaussian elimination */
        double factor;
        double pivot;

        pivot = A->data[k * cols + k];
        if (pivot == 0.0) {
            LOG_ERROR("Error: Zero pivot found at position %zu - Matrix is singular", k);
            return -1;
        }

        for (size_t i = k + 1; i < rows; i++) {
            factor = A->data[i * cols + k] / pivot;
            A->data[i * cols + k] = factor;

            for (size_t j = k + 1; j < cols; j++) {
                A->data[i * cols + j] -= factor * A->data[k * cols + j];
            }
        }
    }

    return 0;
}

int linalg_lu_solve(const matrix_t *LU, const permutation_t *p, const vector_t *v, vector_t *x)
{
    vector_t *y;
    size_t rows;
    size_t cols;

    rows = LU->size1;
    cols = LU->size2;
    y = vector_alloc(rows);

    /* Forward substitution (Ly = Pv) */
    for (size_t i = 0; i < rows; i++) {
        double sum;
        size_t p_i;

        sum = 0.0;
        p_i = p->data[i];
        for (size_t j = 0; j < i; j++) {
            sum += LU->data[i * cols + j] * y->data[j];
        }
        y->data[i] = v->data[p_i] - sum;
    }

    /* Backward substitution (Ux = y) */
    for (int i = (rows - 1); i >= 0; i--) {
        double sum;

        sum = 0.0;
        for (size_t j = (i + 1); j < cols; j++) {
            sum += LU->data[i * cols + j] * x->data[j];
        }
        x->data[i] = (y->data[i] - sum) / LU->data[i * cols + i];
    }

    vector_free(y);
    return 0;
}

int linalg_sv_decomp(matrix_t *A, matrix_t *V, vector_t *S)
{
    size_t rows;
    size_t cols;
    matrix_t *A_copy;
    matrix_t *B;

    rows = A->size1;
    cols = A->size2;

    A_copy = matrix_alloc(rows, cols);
    B = matrix_alloc(cols, cols);

    if (A_copy == NULL || B == NULL)
        return -1;

    // Copy the original matrix A into A_copy
    for (size_t i = 0; i < rows * cols; i++)
        A_copy->data[i] = A->data[i];

    // Compute B = A^T * A
    product_AtA(A_copy, B);

    // Compute eigenvalues and eigenvectors of B → singular values and V
    jacobi_svd(B, V, S);

    // Compute U = A * V * Σ⁻¹ and store the result back into A
    compute_U_from_AV(A_copy, V, S, A);

    matrix_free(A_copy);
    matrix_free(B);

    return 0;
}
