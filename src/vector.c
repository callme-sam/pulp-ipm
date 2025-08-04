#include "vector.h"

vector_t *vector_alloc(const size_t n)
{
    vector_t *v;

    v = malloc(sizeof(vector_t));
    if (!v)
        return NULL;

    v->gsl = gsl_vector_calloc(n);
    if (!v->gsl) {
        free(v);
        return NULL;
    }

    v->size = v->gsl->size;

    return v;
}

void vector_free(vector_t *v)
{
    if (!v) return;

    if (v->gsl)
        gsl_vector_free(v->gsl);
    free(v);
}

void vector_set(vector_t *v, const size_t i, double x)
{
    gsl_vector_set(v->gsl, i, x);
}

double vector_get(const vector_t *v, const size_t i)
{
    return gsl_vector_get(v->gsl, i);
}

int vector_mul(vector_t *a, const vector_t *b)
{
    return gsl_vector_mul(a->gsl, b->gsl);
}

int vector_sub(vector_t *a, const vector_t *b)
{
    return gsl_vector_sub(a->gsl, b->gsl);
}

int vector_add(vector_t *a, const vector_t *b)
{
    return gsl_vector_add(a->gsl, b->gsl);
}

int vector_scale(vector_t *a, const double x)
{
    return gsl_vector_scale(a->gsl, x);
}

int vector_memcpy(vector_t *dst, const vector_t *src)
{
    return gsl_vector_memcpy(dst->gsl, src->gsl);
}

double vector_min(const vector_t *v)
{
    return gsl_vector_min(v->gsl);
}

void vector_set_all(vector_t *v, double x)
{
    return gsl_vector_set_all(v->gsl, x);
}

void vector_set_zero(vector_t *v)
{
    return gsl_vector_set_zero(v->gsl);
}
