#include "log.h"

#include "vector.h"

vector_t *vector_alloc(const size_t n)
{
    vector_t *v;

    v = malloc(sizeof(vector_t));
    if (!v)
        return NULL;

    // TODO: remove - alloc v->data directly
    v->gsl = gsl_vector_calloc(n);
    if (!v->gsl) {
        free(v);
        return NULL;
    }

    v->size = v->gsl->size;
    v->data = v->gsl->data;

    return v;
}

void vector_free(vector_t *v)
{
    if (!v) return;

    // TODO: remove - free v->data directly
    if (v->gsl)
        gsl_vector_free(v->gsl);
    free(v);
}

void vector_set(vector_t *v, const size_t i, double x)
{
    v->data[i] = x;
}

double vector_get(const vector_t *v, const size_t i)
{
    return v->data[i];
}

int vector_mul(vector_t *a, const vector_t *b)
{
    size_t a_len;
    size_t b_len;

    a_len = a->size;
    b_len = b->size;
    if (a_len != b_len) {
        LOG_ERROR("Error: Vectors must have the same length");
        return -1;
    }

    for (size_t i = 0; i < a_len; i++)
        a->data[i] = a->data[i] * b->data[i];

    return 0;
}

int vector_sub(vector_t *a, const vector_t *b)
{
    size_t a_len;
    size_t b_len;

    a_len = a->size;
    b_len = b->size;
    if (a_len != b_len) {
        LOG_ERROR("Error: Vectors must have the same length");
        return -1;
    }

    for (size_t i = 0; i < a_len; i++)
        a->data[i] = a->data[i] - b->data[i];

    return 0;
}

int vector_add(vector_t *a, const vector_t *b)
{
    size_t a_len;
    size_t b_len;

    a_len = a->size;
    b_len = b->size;
    if (a_len != b_len) {
        LOG_ERROR("Error: Vectors must have the same length");
        return -1;
    }

    for (size_t i = 0; i < a_len; i++)
        a->data[i] = a->data[i] + b->data[i];

    return 0;
}

int vector_scale(vector_t *a, const double x)
{
    size_t a_len;

    a_len = a->size;
    for (size_t i = 0; i < a_len; i++)
        a->data[i] = a->data[i] * x;

    return 0;
}

int vector_memcpy(vector_t *dst, const vector_t *src)
{
    size_t src_len;
    size_t dst_len;

    src_len = src->size;
    dst_len = dst->size;
    if (src_len != dst_len) {
        LOG_ERROR("Error: Vectors must have the same length");
        return -1;
    }

    for (size_t i = 0; i < src_len; i++)
        dst->data[i] = src->data[i];

    return 0;
}

double vector_min(const vector_t *v)
{
    size_t v_len;
    double min;

    min = v->data[0];
    v_len = v->size;
    for (size_t i = 1; i < v_len; i++) {
        double x = v->data[i];
        if (x < min)
            min = x;
    }

    return min;
}

void vector_set_all(vector_t *v, double x)
{
    size_t v_len;

    v_len = v->size;
    for (size_t i = 0; i < v_len; i++)
        v->data[i] = x;
}

void vector_set_zero(vector_t *v)
{
    return vector_set_all(v, 0);
}

void vector_set_ones(vector_t *v)
{
    return vector_set_all(v, 1);
}

vector_t *vector_concat(const vector_t *v, const double x)
{
    vector_t *res;
    size_t v_len;

    v_len = v->size;
    res = vector_alloc(v_len + 1);
    if (!res) {
        LOG_ERROR("Error: failed to allocate memory for new vector");
        return NULL;
    }

    for (size_t i = 0; i < v_len; i++)
        res->data[i] = v->data[i];
    res->data[v_len] = x;

    return res;
}

bool has_nonpositive_elements(const vector_t *v)
{
    return vector_min(v) <= 0;
}
