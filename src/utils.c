// #include <stdio.h>

#include "utils.h"

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

void print_vector(const gsl_vector *v)
{
    size_t v_len;

    v_len = v->size;

    for (size_t i = 0; i < v_len; i++)
        printf("%.2f ", gsl_vector_get(v, i));
    printf("\n");
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