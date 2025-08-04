#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "blas.h"
#include "rng.h"
#include "utils.h"

/**
 * @brief Generate a random feasible linear program (LP).
 *
 * Fills the provided matrices/vectors with data such that the LP:
 *     minimize cᵀx  subject to A x = b, x > 0
 * has a known feasible solution.
 *
 * The matrix A is generated with normally distributed entries, the vector c
 * with uniform values in [0,1), and b is computed as b = A * x₀, where x₀ is a
 * strictly positive vector.
 *
 * @param[out] A Pointer to an allocated GSL matrix (m × n) to be filled.
 * @param[out] b Pointer to an allocated GSL vector (m) to be filled.
 * @param[out] c Pointer to an allocated GSL vector (n) to be filled.
 *
 * @note A, b, and c must be allocated before calling this function.
 */
void generate_lp(matrix_t **A, vector_t **b, vector_t **c)
{
    vector_t *tmp;
    size_t rows;
    size_t cols;

    rows = (*A)->size1;
    cols = (*A)->size2;

    assert(rows <= cols);

    tmp = vector_alloc(cols);

    rng_env_setup();
    rng_t *rng = rng_alloc();
    rng_set(rng, time(NULL));

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double val = rand_normal(rng);
            matrix_set(*A, i, j, val);
        }
    }

    // A[0, :] = rand + 0.1 (boundedness)
    for (size_t j = 0; j < cols; j++) {
        double val = rand_uniform(rng) + 0.1;
        matrix_set(*A, 0, j, val);
    }

    for (size_t j = 0; j < cols; j++) {
        double val = rand_uniform(rng) + 0.01;
        vector_set(tmp, j, val);
    }

    // b = A * tmp
    blas_dgemv(BLAS_NO_TRANSPOSE, 1.0, *A, tmp, 0.0, *b);

    // c ~ U(0,1)
    for (size_t j = 0; j < cols; j++) {
        vector_set(*c, j, rand_uniform(rng));
    }

    vector_free(tmp);
    rng_free(rng);
}

/**
 * @brief Fill a GSL matrix with a deterministic pattern based on indices.
 *
 * Each element M(i,j) is set to (i * j + j).
 *
 * @param[out] M Pointer to a GSL matrix to be filled.
 */
void fill_matrix(matrix_t *M)
{
    size_t rows;
    size_t cols;

    rows = M->size1;
    cols = M->size2;

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double val;

            val = (i * j) + j;
            matrix_set(M, i, j, val);
        }
    }
}

/**
 * @brief Fill a GSL vector with sequential values starting from 0.
 *
 * Each element v[i] is set to i.
 *
 * @param[out] v Pointer to a GSL vector to be filled.
 */
void fill_vector(vector_t *v)
{
    size_t v_len;

    v_len = v->size;

    for (size_t i = 0; i < v_len; i++) {
        double val;

        val = i;
        vector_set(v, i, val);
    }

}

/**
 * @brief Print a GSL matrix to standard output in human-readable format.
 *
 * Each row is printed on a new line with values formatted to 2 decimal places.
 *
 * @param[in] M Pointer to the GSL matrix to print.
 */
void print_matrix(const matrix_t *M)
{
    size_t rows;
    size_t cols;

    rows = M->size1;
    cols = M->size2;

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            printf("%.2f ", matrix_get(M, i, j));
        }
        printf("\n");
    }
}

/**
 * @brief Print a GSL vector to standard output in human-readable format.
 *
 * All elements are printed on a single line, formatted to 2 decimal places.
 *
 * @param[in] v Pointer to the GSL vector to print.
 */
void print_vector(const vector_t *v)
{
    size_t v_len;

    v_len = v->size;

    for (size_t i = 0; i < v_len; i++)
        printf("%.2f ", vector_get(v, i));
    printf("\n");
}

/**
 * @brief Logs the contents of a GSL matrix with a custom label.
 *
 * This function prints the contents of a GSL matrix to the log output,
 * formatted as a human-readable 2D array with a descriptive label.
 * Each row is printed on a separate line, with elements separated by commas:
 * ```
 * Matrix <name> (size = rows x cols):
 *     [a_00, a_01, ..., a_0n]
 *     [a_10, a_11, ..., a_1n]
 *     ...
 * ```
 * Logging only occurs if the specified `level` is less than or equal to `CURRENT_LOG_LEVEL`.
 *
 * @param[in] level Logging severity level (e.g., DEBUG, INFO, ERROR).
 * @param[in] m Pointer to the GSL matrix to be logged.
 * @param[in] name Descriptive name of the matrix, shown in the log output.
 *
 * @note The function dynamically allocates memory to format the matrix string.
 *       It automatically releases the memory before returning.
 *
 * @warning If memory allocation fails, the function logs an error and returns without printing.
 */
void log_matrix(LogLevel level, const matrix_t *m, const char *name);

void log_matrix(LogLevel level, const matrix_t *m, const char *name) {
    if (level > CURRENT_LOG_LEVEL) return;

    const size_t per_elem_len = 24;
    const size_t header_len = 64;
    const size_t indent_len = 4;
    size_t row_len = indent_len + per_elem_len * m->size2 + 2;
    size_t buf_size = header_len + row_len * m->size1 + 1;

    char *buf = malloc(buf_size);
    if (!buf) {
        LOG_ERROR("Failed to allocate memory for matrix logging");
        return;
    }

    size_t offset = snprintf(buf, buf_size, "Matrix %s (size = %zux%zu):\n", name, m->size1, m->size2);

    for (size_t i = 0; i < m->size1; ++i) {
        offset += snprintf(buf + offset, buf_size - offset, "\t[");
        for (size_t j = 0; j < m->size2; ++j) {
            offset += snprintf(buf + offset, buf_size - offset, "%g%s",
                               matrix_get(m, i, j), (j < m->size2 - 1) ? ", " : "");
        }
        offset += snprintf(buf + offset, buf_size - offset, "]\n");
    }
    buf[offset] = '\0';

    LOG(level, "%s", buf);
    free(buf);
}

/**
 * @brief Logs the contents of a GSL vector with a custom label.
 *
 * This function formats and prints the elements of a GSL vector using a specified log level.
 * The vector is printed in the format:
 * ```
 * Vector <name> (size = n):
 *     [e0, e1, ..., en]
 * ```
 * Logging only occurs if the specified level is less than or equal to the current log level.
 *
 * @param[in] level Logging severity level
 * @param[in] v Pointer to the GSL vector to be printed.
 * @param[in] name Descriptive name to be shown in the log output.
 *
 * @note The function allocates temporary memory to construct the formatted log message.
 *       It automatically frees the memory before returning.
 *
 * @warning If memory allocation fails, the function logs an error and exits early.
 */
void log_vector(LogLevel level, const vector_t *v, const char *name) {
   if (level > CURRENT_LOG_LEVEL) return;

   const size_t per_elem_len = 24;
    const size_t header_len = 64;
    size_t buf_size = header_len + per_elem_len * v->size;
    char *buf = malloc(buf_size);
    if (!buf) {
        LOG_ERROR("Failed to allocate memory for vector logging");
        return;
    }

    size_t offset = snprintf(buf, buf_size, "Vector %s (size = %zu):\n\t[", name, v->size);

    for (size_t i = 0; i < v->size; ++i) {
        offset += snprintf(buf + offset, buf_size - offset, "%g%s",
                           vector_get(v, i), (i < v->size - 1) ? ", " : "");
    }

    snprintf(buf + offset, buf_size - offset, "]");

    LOG(level, "%s", buf);

    free(buf);
}

/**
 * @brief Checks if the minimum value in a GSL vector is <= 0.
 *
 * @param[in] x Input vector (must not be NULL).
 *
 *
 * @return true if any element in x is <= 0, false otherwise.
 */
bool has_nonpositive_elements(const vector_t *v)
{
    return vector_min(v) <= 0.0;
}

/**
 * @brief Creates a new GSL vector of a specified size and initializes all its elements to 1.0.
 *
 * @param n The size of the vector to be created.
 *
 * @return vector_t* A pointer to the newly allocated and initialized GSL vector,
 * or NULL if memory allocation fails.
 *
 * @note It is the programmer's responsibility to free the allocated
 * memory using `vector_free()` when the vector is no longer needed.
 */
vector_t *vector_ones(size_t n)
{
    vector_t *ones = vector_alloc(n);
    vector_set_all(ones, 1.0);
    return ones;
}

 /**
 * @brief Concatenates a scalar value to the end of a GSL vector.
 *
 * @param v A pointer to the input GSL vector. This vector is not modified.
 * @param b The double-precision scalar value to append to the vector.
 *
 * @return vector_t* A pointer to the newly allocated GSL vector containing the
 * concatenated elements. Returns NULL if memory allocation fails.
 *
 * @note is the programmer's responsibility to free the allocated memory using
 * `vector_free()` when the returned vector is no longer needed.
 */
vector_t *vector_concat(const vector_t *v, double b)
{
    vector_t *res;
    size_t v_len;

    v_len = v->size;
    res = vector_alloc(v_len + 1);

    if (res == NULL)
        return NULL;

    for (size_t i = 0; i < v_len; i++)
        vector_set(res, i, vector_get(v, i));
    vector_set(res, v_len, b);

    return res;
}

const char* err_to_str(const int err)
{
    return gsl_strerror(err);
}