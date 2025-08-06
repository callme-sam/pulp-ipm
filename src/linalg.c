#include <gsl/gsl_linalg.h>

#include "log.h"

#include "linalg.h"

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

int linalg_sv_decomp(matrix_t *A, matrix_t *V, vector_t *S, vector_t *work)
{
    return gsl_linalg_SV_decomp(A->gsl, V->gsl, S->gsl, work->gsl);
}