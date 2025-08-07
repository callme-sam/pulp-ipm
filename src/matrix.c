#include "log.h"

#include "matrix.h"

/**
 * @brief Allocates a new matrix with given dimensions.
 *
 * Allocates memory for a matrix structure and its data array
 * with `n1` rows and `n2` columns.
 *
 * @param[in] n1 Number of rows.
 * @param[in] n2 Number of columns.
 *
 * @return matrix_t* Pointer to the allocated matrix, or NULL on failure.
 *
 * @note The caller is responsible for freeing the memory with permutation_free().
 */
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

/**
 * @brief Frees the memory associated with a matrix.
 *
 * Frees the matrix data and the matrix structure itself.
 * Safe to call with NULL pointer.
 *
 * @param[in,out] m Pointer to the matrix to free.
 */
void matrix_free(matrix_t *m)
{
    if (!m) return;

    if (m->data) free(m->data);
    free(m);
}

/**
 * @brief Sets the element at position (i, j) in the matrix.
 *
 * Stores the value `x` at row `i` and column `j` of matrix `m`.
 * No bounds checking is performed.
 *
 * @param[in,out] m Pointer to the matrix.
 * @param[in] i Row index.
 * @param[in] j Column index.
 * @param[in] x Value to set.
 */
void matrix_set(matrix_t *m, const size_t i, const size_t j, const double x)
{
    size_t n_cols;

    n_cols = m->size2;
    m->data[i * n_cols + j] = x;
}

/**
 * @brief Gets the element at position (i, j) from the matrix.
 *
 * Returns the value stored at row `i` and column `j` of matrix `m`.
 * No bounds checking is performed.
 *
 * @param[in] m Pointer to the matrix.
 * @param[in] i Row index.
 * @param[in] j Column index.
 *
 * @return double The element at position (i, j).
 */
double matrix_get(const matrix_t *m, const size_t i, const size_t j)
{
    size_t n_cols;
    double x;

    n_cols = m->size2;
    x = m->data[i * n_cols + j];

    return x;
}

/**
 * @brief Copies the transpose of `src` into `dst`.
 *
 * Copies the transposed contents of matrix `src` into matrix `dst`.
 * The sizes of `dst` and `src` must match the transpose dimensions.
 *
 * @param[out] dst Pointer to the destination matrix (transposed).
 * @param[in] src Pointer to the source matrix.
 *
 * @return int 0 on success, -1 if dimensions do not match.
 */
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

/**
 * @brief Copies the contents of one matrix to another.
 *
 * Copies all elements from `src` to `dst`. The matrices must have the same size.
 *
 * @param[out] dst Pointer to the destination matrix.
 * @param[in] src Pointer to the source matrix.
 *
 * @return int 0 on success, -1 if sizes differ.
 */
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

/**
 * @brief Swaps two rows in a matrix.
 *
 * Exchanges the contents of row `i` and row `j` in matrix `m`.
 *
 * @param[in,out] m Pointer to the matrix.
 * @param[in] i Index of the first row.
 * @param[in] j Index of the second row.
 */
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

/**
 * @brief Sets all elements of a matrix to 1.0.
 *
 * Fills the matrix `m` with ones.
 *
 * @param[in,out] m Pointer to the matrix to set.
 */
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