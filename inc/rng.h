#ifndef RNG_H_
#define RNG_H_

#include <gsl/gsl_rng.h>

typedef struct
{
    unsigned long int state;

    // TODO: remove
    gsl_rng *gsl;
} rng_t;

rng_t *rng_alloc();
void rng_free(rng_t *r);

void rng_set(rng_t *r, unsigned long int seed);
double rng_normal(rng_t *rng);
double rng_uniform(rng_t *r);


#endif  /* RNG_H_ */