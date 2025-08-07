#include "log.h"

#include "blas.h"

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