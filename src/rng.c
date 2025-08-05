#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <gsl/gsl_randist.h>

#include "rng.h"

rng_t *rng_alloc()
{
    rng_t *r;

    r = malloc(sizeof(rng_t));
    if (!r)
        return NULL;

    r->state = 1;

    // TODO: remove
    r->gsl = gsl_rng_alloc(gsl_rng_default);
    if (!r->gsl) {
        free(r);
        return NULL;
    }

    return r;
}

void rng_free(rng_t *r)
{
    if (!r) return;

    // TODO: remove
    if (r->gsl)
        gsl_rng_free(r->gsl);
    free(r);    // keep
}

void rng_set(rng_t *r, unsigned long int seed)
{
    r->state = (seed == 0) ? 1 : seed;
}

double rng_uniform(rng_t *r)
{
    // Uses XorShift @ 32-bit
    uint32_t x;

    x = (uint32_t) r->state;
    x = x ^ (x << 13);
    x = x ^ (x >> 17);
    x = x ^ (x << 5);

    r->state = x;
    return (double) x / (double) UINT32_MAX;
}

double rng_normal(rng_t *rng)
{
    // Uses Box-Muller method
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