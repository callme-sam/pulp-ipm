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