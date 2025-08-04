#include <gsl/gsl_randist.h>

#include "rng.h"

rng_t *rng_alloc()
{
    rng_t *r;

    r = malloc(sizeof(rng_t));
    if (!r)
        return NULL;

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

    if (r->gsl)
        gsl_rng_free(r->gsl);
    free(r);
}

void rng_set(const rng_t *r, unsigned long int seed)
{
    gsl_rng_set(r->gsl, seed);
}

double rng_uniform(const rng_t *r)
{
    return gsl_rng_uniform(r->gsl);
}

void rng_env_setup()
{
    gsl_rng_env_setup();
}

double rand_uniform(rng_t *rng)
{
    return gsl_rng_uniform(rng->gsl);
}

double rand_normal(rng_t *rng)
{
    return gsl_ran_gaussian(rng->gsl, 1.0);
}