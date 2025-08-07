#include "permutation.h"

/**
 * @brief Allocates memory for a new permutation vector.
 *
 * Creates a permutation vector of size `n` initialized with undefined values.
 * Memory is allocated for both the permutation structure and its data array.
 *
 * @param[in] n Size of the permutation vector.
 * @return permutation_t* Pointer to the newly allocated permutation vector.
 * @retval NULL if memory allocation fails.
 *
 * @note The caller is responsible for freeing the memory with permutation_free().
 */
permutation_t *permutation_alloc(const size_t n)
{
    permutation_t *p;

    p = malloc(sizeof(permutation_t));
    if (!p)
        return NULL;

    p->size = n;
    p->data = malloc(n * sizeof(double));
    if (!p->data) {
        free(p);
        return NULL;
    }

    return p;
}

/**
 * @brief Frees memory allocated for a permutation vector.
 *
 * Safely deallocates the memory of both the permutation data array and the structure itself.
 * Handles NULL pointers gracefully.
 *
 * @param[in,out] p Pointer to the permutation vector to free.
 *
 * @note This function is safe to call with NULL pointers.
 */
void permutation_free(permutation_t *p)
{
    if (!p) return;

    if (p->data) free(p->data);
    free(p);
}

/**
 * @brief Initializes a permutation vector as the identity permutation.
 *
 * Sets each element p[i] = i, creating a no-op permutation.
 *
 * @param[in,out] p Permutation vector to initialize.
 *
 * @pre p must be a valid permutation vector (non-NULL and properly allocated).
 * @post After execution, p represents the identity permutation.
 */
void permutation_set_identity(permutation_t *p)
{
    size_t p_len;

    p_len = p->size;
    for (size_t i = 0; i < p_len; i++)
        p->data[i] = i;
}

/**
 * @brief Swaps two elements in a permutation vector.
 *
 * Exchanges the values at positions i and j in the permutation vector.
 * Used primarily for row swapping in matrix factorizations.
 *
 * @param[in,out] p Permutation vector to modify.
 * @param[in] i First index to swap (0-based).
 * @param[in] j Second index to swap (0-based).
 *
 * @pre Both i and j must be valid indices (0 <= i,j < p->size).
 * @pre p must be non-NULL and properly allocated.
 */
void permutation_swap(permutation_t *p, const size_t i, const size_t j)
{
    size_t tmp;

    tmp = p->data[i];
    p->data[i] = p->data[j];
    p->data[j] = tmp;
}