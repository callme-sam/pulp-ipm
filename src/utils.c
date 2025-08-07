#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "blas.h"
#include "rng.h"

#include "utils.h"

/**
 * @brief Generates a random feasible linear programming (LP) problem.
 *
 * Constructs a well-conditioned LP problem in standard form:
 * \[
 * \begin{aligned}
 * \text{minimize} \quad & c^T x \\
 * \text{subject to} \quad & A x = b \\
 * & x \geq 0
 * \end{aligned}
 * \]
 * with known feasible solution x₀. The problem is generated such that:
 * - Matrix A has normally distributed entries (except first row)
 * - First row of A ensures boundedness (uniform entries + 0.1)
 * - Vector b is computed as A * x₀ where x₀ has uniform entries in [0.01, 1.01)
 * - Cost vector c has uniform entries in [0,1)
 *
 * @param[in,out] A Pointer to pre-allocated matrix (m×n) to fill with constraint coefficients.
 *                   Will be modified to contain the generated matrix.
 * @param[in,out] b Pointer to pre-allocated vector (m) to fill with RHS constraints.
 *                   Will be modified to contain b = A * x₀.
 * @param[in,out] c Pointer to pre-allocated vector (n) to fill with cost coefficients.
 *                   Will be modified with uniform random costs.
 *
 * @pre Matrices/Vectors must be pre-allocated with correct dimensions
 * @pre rows <= cols (system must not be overconstrained)
 * @post A, b, c contain a feasible LP problem with known solution x₀ = tmp
 * @post The generated problem satisfies strict feasibility (x₀ > 0)
 *
 * @note Uses XorShift for uniform numbers and Box-Muller for normal distribution
 * @warning Seed is based on current time - not suitable for cryptographic purposes
 * @warning Assumes matrices/vectors are properly allocated (will assert on rows > cols)
 *
 * @example
 * // Generate a 10x20 LP problem
 * matrix_t *A = matrix_alloc(10, 20);
 * vector_t *b = vector_alloc(10);
 * vector_t *c = vector_alloc(20);
 * generate_lp(&A, &b, &c);
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

    rng_t *rng = rng_alloc();
    rng_set(rng, time(NULL));

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double val = rng_normal(rng);
            matrix_set(*A, i, j, val);
        }
    }

    // A[0, :] = rand + 0.1 (boundedness)
    for (size_t j = 0; j < cols; j++) {
        double val = rng_uniform(rng) + 0.1;
        matrix_set(*A, 0, j, val);
    }

    for (size_t j = 0; j < cols; j++) {
        double val = rng_uniform(rng) + 0.01;
        vector_set(tmp, j, val);
    }

    // b = A * tmp
    blas_dgemv(BLAS_NO_TRANSPOSE, 1.0, *A, tmp, 0.0, *b);

    // c ~ U(0,1)
    for (size_t j = 0; j < cols; j++) {
        vector_set(*c, j, rng_uniform(rng));
    }

    vector_free(tmp);
    rng_free(rng);
}