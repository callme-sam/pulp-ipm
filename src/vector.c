#include "log.h"

#include "vector.h"

/**
 * @brief Allocates a new vector of floating-point numbers.
 *
 * This function dynamically allocates memory for a `vector_t` structure,
 * which represents a vector of `n` double-precision floating-point elements.
 * It first allocates the structure itself, then allocates memory for the
 * data array. If either allocation fails, the function returns `NULL`.
 *
 * @param[in] n The number of elements to allocate in the vector.
 *
 * @return vector_t* A pointer to the newly allocated vector, or NULL on failure.
 *
 * @note The caller is responsible for freeing the memory with permutation_free().
 */
vector_t *vector_alloc(const size_t n)
{
    vector_t *v;

    v = malloc(sizeof(vector_t));
    if (!v)
        return NULL;

    v->size = n;
    v->data = malloc(n * sizeof(double));
    if (!v->data) {
        free(v);
        return NULL;
    }

    return v;
}

/**
 * @brief Frees the memory allocated for a vector.
 *
 * This function deallocates the memory used by a `vector_t` structure,
 * including its internal data array. If the input pointer is NULL,
 * or if the data pointer inside the structure is NULL, the function
 * safely does nothing.
 *
 * @param[in] v Pointer to the vector to free. Can be NULL.
 */
void vector_free(vector_t *v)
{
    if (!v) return;

    if (v->data) free(v->data);
    free(v);
}

/**
 * @brief Sets the value of an element in a vector.
 *
 * This function assigns the value `x` to the `i`-th element of the vector `v`.
 * It does not perform bounds checking, so the caller must ensure that `i` is
 * less than `v->size` and that `v` and `v->data` are not NULL.
 *
 * @param[in,out] v Pointer to the vector. The function modifies its internal data.
 * @param[in]     i Index of the element to set (must be in range [0, v->size - 1]).
 * @param[in]     x The value to assign to the element at index `i`.
 */
void vector_set(vector_t *v, const size_t i, double x)
{
    v->data[i] = x;
}

/**
 * @brief Retrieves the value of an element in a vector.
 *
 * This function returns the value stored at index `i` in the vector `v`.
 * No bounds checking is performed, so the caller must ensure that `i` is
 * within the valid range and that `v` and `v->data` are not NULL.
 *
 * @param[in] v Pointer to the vector from which to retrieve the value.
 * @param[in] i Index of the element to retrieve (must be in range [0, v->size - 1]).
 *
 * @return double The value stored at index `i` in the vector.
 */
double vector_get(const vector_t *v, const size_t i)
{
    return v->data[i];
}

/**
 * @brief Performs element-wise multiplication of two vectors.
 *
 * This function multiplies each element of vector `a` by the corresponding element
 * of vector `b`, storing the result back in vector `a`. Both vectors must have the
 * same size. If the sizes do not match, the function logs an error and returns -1.
 *
 * @param[in,out] a Pointer to the first vector. It will be modified in-place to store the result.
 * @param[in]     b Pointer to the second vector. It is read-only and must be the same length as `a`.
 *
 * @return int Returns 0 on success, or -1 if the vectors have mismatched sizes.
 */
int vector_mul(vector_t *a, const vector_t *b)
{
    size_t a_len;
    size_t b_len;

    a_len = a->size;
    b_len = b->size;
    if (a_len != b_len) {
        LOG_ERROR("Error: Vectors must have the same length");
        return -1;
    }

    for (size_t i = 0; i < a_len; i++)
        a->data[i] = a->data[i] * b->data[i];

    return 0;
}

/**
 * @brief Performs element-wise subtraction of two vectors.
 *
 * Subtracts each element of vector `b` from the corresponding element in vector `a`,
 * and stores the result in vector `a`. Both vectors must have the same length.
 *
 * @param[in,out] a Pointer to the vector that will store the result.
 * @param[in]     b Pointer to the vector to subtract from `a`.
 *
 * @return int Returns 0 on success, -1 if vector lengths do not match.
 */
int vector_sub(vector_t *a, const vector_t *b)
{
    size_t a_len;
    size_t b_len;

    a_len = a->size;
    b_len = b->size;
    if (a_len != b_len) {
        LOG_ERROR("Error: Vectors must have the same length");
        return -1;
    }

    for (size_t i = 0; i < a_len; i++)
        a->data[i] = a->data[i] - b->data[i];

    return 0;
}

/**
 * @brief Performs element-wise addition of two vectors.
 *
 * Adds each element of vector `b` to the corresponding element in vector `a`,
 * storing the result in vector `a`. Both vectors must have the same length.
 *
 * @param[in,out] a Pointer to the vector that will store the result.
 * @param[in]     b Pointer to the vector to add to `a`.
 *
 * @return int Returns 0 on success, -1 if vector lengths do not match.
 */
int vector_add(vector_t *a, const vector_t *b)
{
    size_t a_len;
    size_t b_len;

    a_len = a->size;
    b_len = b->size;
    if (a_len != b_len) {
        LOG_ERROR("Error: Vectors must have the same length");
        return -1;
    }

    for (size_t i = 0; i < a_len; i++)
        a->data[i] = a->data[i] + b->data[i];

    return 0;
}

/**
 * @brief Scales all elements of a vector by a constant factor.
 *
 * Multiplies each element of vector `a` by the scalar `x`, modifying `a` in-place.
 *
 * @param[in,out] a Pointer to the vector to scale.
 * @param[in]     x Scalar multiplier.
 *
 * @return int Always returns 0.
 */
int vector_scale(vector_t *a, const double x)
{
    size_t a_len;

    a_len = a->size;
    for (size_t i = 0; i < a_len; i++)
        a->data[i] = a->data[i] * x;

    return 0;
}

/**
 * @brief Adds a constant to all elements of a vector.
 *
 * Adds the scalar value `x` to each element of vector `a`, modifying `a` in-place.
 *
 * @param[in,out] a Pointer to the vector to modify.
 * @param[in]     x Scalar value to add.
 *
 * @return int Always returns 0.
 */
int vector_add_constant(vector_t *a, const double x)
{
    size_t a_len;

    a_len = a->size;
    for (size_t i = 0; i < a_len; i++)
        a->data[i] = a->data[i] + x;

    return 0;
}

/**
 * @brief Copies the contents of one vector into another.
 *
 * Copies all elements from vector `src` to vector `dst`. Both vectors must have the same size.
 *
 * @param[out]    dst Pointer to the destination vector.
 * @param[in]     src Pointer to the source vector.
 *
 * @return int Returns 0 on success, -1 if vector sizes do not match.
 */
int vector_memcpy(vector_t *dst, const vector_t *src)
{
    size_t src_len;
    size_t dst_len;

    src_len = src->size;
    dst_len = dst->size;
    if (src_len != dst_len) {
        LOG_ERROR("Error: Vectors must have the same length");
        return -1;
    }

    for (size_t i = 0; i < src_len; i++)
        dst->data[i] = src->data[i];

    return 0;
}

/**
 * @brief Returns the minimum value in a vector.
 *
 * Scans the vector `v` and returns the smallest element.
 *
 * @param[in] v Pointer to the vector.
 *
 * @return double The minimum value in the vector.
 */
double vector_min(const vector_t *v)
{
    size_t v_len;
    double min;

    min = v->data[0];
    v_len = v->size;
    for (size_t i = 1; i < v_len; i++) {
        double x = v->data[i];
        if (x < min)
            min = x;
    }

    return min;
}

/**
 * @brief Sets all elements of a vector to a given value.
 *
 * Assigns the value `x` to every element of vector `v`.
 *
 * @param[in,out] v Pointer to the vector to modify.
 * @param[in]     x Value to assign to each element.
 */
void vector_set_all(vector_t *v, double x)
{
    size_t v_len;

    v_len = v->size;
    for (size_t i = 0; i < v_len; i++)
        v->data[i] = x;
}

/**
 * @brief Sets all elements of a vector to zero.
 *
 * Shortcut for calling `vector_set_all(v, 0)`.
 *
 * @param[in,out] v Pointer to the vector to zero.
 */
void vector_set_zero(vector_t *v)
{
    return vector_set_all(v, 0);
}

/**
 * @brief Sets all elements of a vector to one.
 *
 * Shortcut for calling `vector_set_all(v, 1)`.
 *
 * @param[in,out] v Pointer to the vector to fill with ones.
 */
void vector_set_ones(vector_t *v)
{
    return vector_set_all(v, 1);
}

/**
 * @brief Returns a new vector by appending a value to an existing vector.
 *
 * Allocates and returns a new vector that contains all elements of vector `v`
 * followed by the value `x`. The original vector `v` is not modified.
 *
 * @param[in] v Pointer to the original vector.
 * @param[in] x Value to append at the end of the new vector.
 *
 * @return vector_t* Pointer to the newly allocated vector, or NULL on failure.
 */
vector_t *vector_concat(const vector_t *v, const double x)
{
    vector_t *res;
    size_t v_len;

    v_len = v->size;
    res = vector_alloc(v_len + 1);
    if (!res) {
        LOG_ERROR("Error: failed to allocate memory for new vector");
        return NULL;
    }

    for (size_t i = 0; i < v_len; i++)
        res->data[i] = v->data[i];
    res->data[v_len] = x;

    return res;
}

/**
 * @brief Checks whether a vector contains any non-positive elements.
 *
 * Returns true if any element in the vector `v` is less than or equal to zero.
 *
 * @param[in] v Pointer to the vector to check.
 *
 * @return true If the vector contains at least one non-positive element.
 * @return false If all elements are strictly positive.
 */
bool has_nonpositive_elements(const vector_t *v)
{
    return vector_min(v) <= 0;
}
