#ifndef PERMUTATION_H_
#define PERMUTATION_H_

#include <gsl/gsl_permutation.h>

typedef struct
{
    size_t size;
    size_t *data;

    // TODO: remove
    gsl_permutation *gsl;
} permutation_t;

permutation_t *permutation_alloc(const size_t n);
void permutation_free(permutation_t *p);

#endif  /* PERMUTATION_H_ */