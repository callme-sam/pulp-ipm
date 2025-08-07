#ifndef PERMUTATION_H_
#define PERMUTATION_H_

#include <stdlib.h>

typedef struct
{
    size_t size;
    size_t *data;
} permutation_t;

permutation_t *permutation_alloc(const size_t n);
void permutation_free(permutation_t *p);

void permutation_set_identity(permutation_t *p);
void permutation_swap(permutation_t *p, const size_t i, const size_t j);

#endif  /* PERMUTATION_H_ */