#include "permutation.h"

permutation_t *permutation_alloc(const size_t n)
{
    permutation_t *p;

    p = malloc(sizeof(permutation_t));
    if (!p)
        return NULL;

    // TODO: remove - alloc p->data directrly
    p->gsl = gsl_permutation_alloc(n);
    if (!p->gsl) {
        free(p);
        return NULL;
    }

    p->size = p->gsl->size;
    p->data = p->gsl->data;

    return p;
}

void permutation_free(permutation_t *p)
{
    if (!p) return;

    // TODO: remove - free p->data directly
    if (p->gsl)
        gsl_permutation_free(p->gsl);
    free(p);
}