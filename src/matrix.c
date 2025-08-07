#include "log.h"

#include "matrix.h"

matrix_t *matrix_alloc(const size_t n1, const size_t n2)
{
    matrix_t *m;

    m = malloc(sizeof(matrix_t));
    if (!m)
        return NULL;

    m->size1 = n1;
    m->size2 = n2;
    m->data = malloc(n1 * n2 * sizeof(double));
    if (!m->data) {
        free(m);
        return NULL;
    }

    return m;
}

void matrix_free(matrix_t *m)
{
    if (!m) return;

    if (m->data) free(m->data);
    free(m);
}

void matrix_set(matrix_t *m, const size_t i, const size_t j, const double x)
{
    size_t n_cols;

    n_cols = m->size2;
    m->data[i * n_cols + j] = x;
}

double matrix_get(const matrix_t *m, const size_t i, const size_t j)
{
    size_t n_cols;
    double x;

    n_cols = m->size2;
    x = m->data[i * n_cols + j];

    return x;
}

int matrix_transpose_memcpy(matrix_t *dst, const matrix_t *src)
{
    size_t src_rows;
    size_t src_cols;
    size_t dst_rows;
    size_t dst_cols;

    src_rows = src->size1;
    src_cols = src->size2;
    dst_rows = dst->size1;
    dst_cols = dst->size2;

    if (src_rows != dst_cols || src_cols != dst_rows) {
        LOG_ERROR("Error: matrix sizes are different");
        return -1;
    }

    for (size_t i = 0; i < src_rows; i++) {
        for (size_t j = 0; j < src_cols; j++) {
            dst->data[i * dst_cols + j] = src->data[j * src_cols + i];
        }
    }

    return 0;
}

int matrix_memcpy(matrix_t *dst, const matrix_t *src)
{
    size_t src_rows;
    size_t src_cols;
    size_t dst_rows;
    size_t dst_cols;

    src_rows = src->size1;
    src_cols = src->size2;
    dst_rows = dst->size1;
    dst_cols = dst->size2;

    if (src_rows != dst_rows || src_cols != dst_cols) {
        LOG_ERROR("Error: matrix sizes are different");
        return -1;
    }

    for (size_t i = 0; i < src_rows; i++) {
        for (size_t j = 0; j < src_cols; j++) {
            dst->data[i * dst_cols + j] = src->data[i * src_cols + j];
        }
    }

    return 0;
}

void matrix_swap_rows(matrix_t *m, const size_t i, const size_t j)
{
    double *row_i;
    double *row_j;
    size_t cols;

    cols = m->size2;
    row_i = &(m->data[i * cols]);
    row_j = &(m->data[j * cols]);

    for (size_t k = 0; k < cols; k++) {
        double tmp = row_i[k];
        row_i[k] = row_j[k];
        row_j[k] = tmp;
    }
}

void matrix_set_identity(matrix_t *m)
{
    size_t rows;
    size_t cols;

    rows = m->size1;
    cols = m->size2;

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            m->data[i * cols + j] = 1.0;
        }
    }
}