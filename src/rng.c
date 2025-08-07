#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "rng.h"

/**
 * @brief Allocates and initializes a new RNG instance.
 *
 * Creates a new random number generator with default initial state (seed = 1).
 *
 * @return rng_t* Pointer to newly allocated RNG instance
 * @retval NULL if memory allocation fails
 *
 * @note The generator should be seeded with rng_set() for proper randomization.
 * @warning The caller is responsible for freeing the RNG with rng_free().
 */
rng_t *rng_alloc()
{
    rng_t *r;

    r = malloc(sizeof(rng_t));
    if (!r)
        return NULL;

    r->state = 1;

    return r;
}

/**
 * @brief Releases memory allocated for an RNG instance.
 *
 * Safely deallocates the RNG object. Handles NULL pointers gracefully.
 *
 * @param[in] r RNG instance to free
 */
void rng_free(rng_t *r)
{
    if (!r) return;
    free(r);
}

/**
 * @brief Seeds the random number generator.
 *
 * Initializes the RNG state. A seed of 0 is automatically converted to 1.
 *
 * @param[in,out] r RNG instance to seed
 * @param[in] seed Seed value (0 is converted to 1)
 *
 * @pre r must be a valid RNG instance (non-NULL)
 * @post The RNG state is initialized for pseudo-random sequence generation
 */
void rng_set(rng_t *r, unsigned long int seed)
{
    r->state = (seed == 0) ? 1 : seed;
}

/**
 * @brief Generates a uniform random number in [0,1).
 *
 * Uses 32-bit XorShift algorithm for high-performance uniform random number generation.
 *
 * @param[in] r Initialized RNG instance
 * @return double Random number in the interval [0,1)
 *
 * @pre r must be properly initialized (via rng_alloc() and optionally rng_set())
 * @note The sequence period is 2^32 - 1
 * @warning Not cryptographically secure
 */
double rng_uniform(rng_t *r)
{
    uint32_t x;

    x = (uint32_t) r->state;
    x = x ^ (x << 13);
    x = x ^ (x >> 17);
    x = x ^ (x << 5);

    r->state = x;
    return (double) x / (double) UINT32_MAX;
}

/**
 * @brief Generates a normally distributed random number (mean=0, variance=1).
 *
 * Uses the Box-Muller transform to convert uniform random numbers to normal distribution.
 *
 * @param[in] rng Initialized RNG instance
 * @return double Random number from standard normal distribution N(0,1)
 *
 * @pre rng must be properly initialized
 * @note Uses polar form of Box-Muller for numerical stability
 */
double rng_normal(rng_t *rng)
{
    double u1;
    double u2;
    double r;
    double f;
    double n;

    do {
        u1 = 2.0 * rng_uniform(rng) - 1.0;  // [-1; 1]
        u2 = 2.0 * rng_uniform(rng) - 1.0;
        r = (u1 * u1) + (u2 * u2);
    } while (r == 0.0 || r >= 1.0);

    f = sqrt(-2.0 * log(r) / r);
    n = u1 * f;

    return n;
}