#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

#include "utils.h"

/**
 * @brief Generate a random number uniformly in the interval [0, 1).
 *
 * @param[in] rng Pointer to a GSL random number generator.
 * @return A random double sampled from the uniform distribution in [0, 1).
 */
static double rand_uniform(gsl_rng *rng)
{
    return gsl_rng_uniform(rng);
}

/**
 * @brief Generate a random number from the standard normal distribution.
 *
 * @param[in] rng Pointer to a GSL random number generator.
 * @return A random double sampled from N(0, 1).
 */
static double rand_normal(gsl_rng *rng) {
    return gsl_ran_gaussian(rng, 1.0);
}

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
void generate_lp(gsl_matrix **A, gsl_vector **b, gsl_vector **c)
{
    gsl_vector *tmp;
    size_t rows;
    size_t cols;

    rows = (*A)->size1;
    cols = (*A)->size2;

    assert(rows <= cols);

    tmp = gsl_vector_alloc(cols);

    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, time(NULL));

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double val = rand_normal(rng);
            gsl_matrix_set(*A, i, j, val);
        }
    }

    // A[0, :] = rand + 0.1 (boundedness)
    for (size_t j = 0; j < cols; j++) {
        double val = rand_uniform(rng) + 0.1;
        gsl_matrix_set(*A, 0, j, val);
    }

    for (size_t j = 0; j < cols; j++) {
        double val = rand_uniform(rng) + 0.01;
        gsl_vector_set(tmp, j, val);
    }

    // b = A * tmp
    gsl_blas_dgemv(CblasNoTrans, 1.0, *A, tmp, 0.0, *b);

    // c ~ U(0,1)
    for (size_t j = 0; j < cols; j++) {
        gsl_vector_set(*c, j, rand_uniform(rng));
    }

    gsl_vector_free(tmp);
    gsl_rng_free(rng);
}

/**
 * @brief Fill a GSL matrix with a deterministic pattern based on indices.
 *
 * Each element M(i,j) is set to (i * j + j).
 *
 * @param[out] M Pointer to a GSL matrix to be filled.
 */
void fill_matrix(gsl_matrix *M)
{
    size_t rows;
    size_t cols;

    rows = M->size1;
    cols = M->size2;

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double val;

            val = (i * j) + j;
            gsl_matrix_set(M, i, j, val);
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
void fill_vector(gsl_vector *v)
{
    size_t v_len;

    v_len = v->size;

    for (size_t i = 0; i < v_len; i++) {
        double val;

        val = i;
        gsl_vector_set(v, i, val);
    }

}

/**
 * @brief Print a GSL matrix to standard output in human-readable format.
 *
 * Each row is printed on a new line with values formatted to 2 decimal places.
 *
 * @param[in] M Pointer to the GSL matrix to print.
 */
void print_matrix(const gsl_matrix *M)
{
    size_t rows;
    size_t cols;

    rows = M->size1;
    cols = M->size2;

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            printf("%.2f ", gsl_matrix_get(M, i, j));
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
void print_vector(const gsl_vector *v)
{
    size_t v_len;

    v_len = v->size;

    for (size_t i = 0; i < v_len; i++)
        printf("%.2f ", gsl_vector_get(v, i));
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
void log_matrix(LogLevel level, const gsl_matrix *m, const char *name);

void log_matrix(LogLevel level, const gsl_matrix *m, const char *name) {
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
                               gsl_matrix_get(m, i, j), (j < m->size2 - 1) ? ", " : "");
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
void log_vector(LogLevel level, const gsl_vector *v, const char *name) {
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
                           gsl_vector_get(v, i), (i < v->size - 1) ? ", " : "");
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
bool has_nonpositive_elements(const gsl_vector *v)
{
    return gsl_vector_min(v) <= 0.0;
}

/**
 * @brief Creates a new GSL vector of a specified size and initializes all its elements to 1.0.
 *
 * @param n The size of the vector to be created.
 *
 * @return gsl_vector* A pointer to the newly allocated and initialized GSL vector,
 * or NULL if memory allocation fails.
 *
 * @note It is the programmer's responsibility to free the allocated
 * memory using `gsl_vector_free()` when the vector is no longer needed.
 */
gsl_vector *vector_ones(size_t n)
{
    gsl_vector *ones = gsl_vector_alloc(n);
    gsl_vector_set_all(ones, 1.0);
    return ones;
}

 /**
 * @brief Concatenates a scalar value to the end of a GSL vector.
 *
 * @param v A pointer to the input GSL vector. This vector is not modified.
 * @param b The double-precision scalar value to append to the vector.
 *
 * @return gsl_vector* A pointer to the newly allocated GSL vector containing the
 * concatenated elements. Returns NULL if memory allocation fails.
 *
 * @note is the programmer's responsibility to free the allocated memory using
 * `gsl_vector_free()` when the returned vector is no longer needed.
 */
gsl_vector *vector_concat(const gsl_vector *v, double b)
{
    gsl_vector *res;
    size_t v_len;

    v_len = v->size;
    res = gsl_vector_alloc(v_len + 1);

    if (res == NULL)
        return NULL;

    for (size_t i = 0; i < v_len; i++)
        gsl_vector_set(res, i, gsl_vector_get(v, i));
    gsl_vector_set(res, v_len, b);

    return res;
}