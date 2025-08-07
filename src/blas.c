#include "log.h"

#include "blas.h"

/**
 * @brief Generalized matrix-vector multiplication: y = α·op(A)·x + β·y
 *
 * Computes one of:
 * - y = α·A·x + β·y   (trans = BLAS_NO_TRANSPOSE)
 * - y = α·Aᵀ·x + β·y  (trans = BLAS_TRANSPOSE)
 *
 * @param[in] trans Whether to transpose A (BLAS_NO_TRANSPOSE/BLAS_TRANSPOSE)
 * @param[in] alpha Scalar multiplier for matrix-vector product
 * @param[in] A Input matrix (size m×n)
 * @param[in] x Input vector (size n if no trans, size m if trans)
 * @param[in] beta Scalar multiplier for y
 * @param[in,out] y Output vector (size m if no trans, size n if trans)
 * @return int 0 on success, -1 on dimension mismatch
 *
 * @note Optimized for alpha=0 and beta=0/1 cases
 * @warning Performs full dimension checking with error logging
 *
 * @example
 * // Compute y = 2.5*A*x + 0.5*y
 * blas_dgemv(BLAS_NO_TRANSPOSE, 2.5, A, x, 0.5, y);
 */
int blas_dgemv(blas_trans_e trans, double alpha, const matrix_t *A, const vector_t *x, double beta, vector_t *y)
{
    size_t rows;
    size_t cols;

    rows = A->size1;
    cols = A->size2;

    /* Optimization */
    if (alpha == 0.0) {
        if (beta == 0.0) {
            vector_set_zero(y);
        } else if (beta != 1.0) {
            vector_scale(y, beta);
        }

        return 0;
    }

    if (trans == BLAS_NO_TRANSPOSE) {
        if (y->size != rows || x->size != cols) {
            LOG_ERROR("Error: dimension mismatch");
            return -1;
        }

        /* y = beta * y + alpha * A * x */
        for (size_t i = 0; i < rows; i++) {
            double sum = 0;
            for (size_t j = 0; j < cols; j++) {
                sum += A->data[i * cols + j] * x->data[j];
            }
            y->data[i] = beta * y->data[i] + alpha * sum;
        }
    } else {
        if (y->size != cols || x->size != rows) {
            LOG_ERROR("Error: dimension mismatch");
            return -1;
        }

        /* y = beta * y + alpha * A^T * x */
        for (size_t i = 0; i < cols; i++) {
            double sum = 0;
            for (size_t j = 0; j < rows; j++) {
                sum += A->data[j * cols +i] * x->data[j];
            }
            y->data[i] = alpha * sum + beta * y->data[i];
        }
    }

    return 0;
}


/**
 * @brief Computes the dot product of two vectors: result = xᵀy
 *
 * Calculates the inner product ∑(x[i]·y[i]) with full error checking.
 *
 * @param[in] x First input vector
 * @param[in] y Second input vector
 * @param[out] result Pointer to store the dot product result
 * @return int 0 on success, -1 if vectors have different lengths
 *
 * @note Uses straightforward accumulation for numerical stability
 * @warning Checks vector lengths and logs errors on mismatch
 */
int blas_ddot(const vector_t *x, const vector_t *y, double *result)
{
    size_t x_len;
    size_t y_len;

    x_len = x->size;
    y_len = y->size;
    if (x_len != y_len) {
        LOG_ERROR("Error: Vectors must have the same length");
        return -1;
    }

    *result = 0;
    for (size_t i = 0; i < x_len; i++)
        *result += x->data[i] * y->data[i];

    return 0;
}

/**
 * @brief Scaled vector addition: y = α·x + y
 *
 * Performs the operation y[i] += α·x[i] for all elements.
 * Includes optimizations for α=0 and full dimension checking.
 *
 * @param[in] alpha Scalar multiplier
 * @param[in] x Input vector to scale
 * @param[in,out] y Vector to add to (modified in-place)
 * @return int 0 on success, -1 if vectors have different lengths
 *
 * @note Skips computation when alpha=0
 * @warning Verifies vector lengths match before operation
 */
int blas_daxpy(double alpha, const vector_t *x, vector_t *y)
{
    size_t x_len;
    size_t y_len;

    x_len = x->size;
    y_len = y->size;
    if (x_len != y_len) {
        LOG_ERROR("Error: Vectors must have the same length");
        return -1;
    }

    if (alpha == 0.0)
        return 0;

    for (size_t i = 0; i < x_len; i++)
        y->data[i] += alpha * x->data[i];

    return 0;
}