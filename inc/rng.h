#ifndef RNG_H_
#define RNG_H_

#include <gsl/gsl_rng.h>

typedef struct
{
    gsl_rng *gsl;
} rng_t;

rng_t *rng_alloc();
void rng_free(rng_t *r);

void rng_set(const rng_t *r, unsigned long int seed);
double rng_uniform(const rng_t *r);
void rng_env_setup();

double rand_uniform(rng_t *rng);
double rand_normal(rng_t *rng);

#endif  /* RNG_H_ */