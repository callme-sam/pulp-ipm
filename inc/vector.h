#ifndef VECTOR_H_
#define VECTOR_H_

#include <stdbool.h>
#include <stdlib.h>

#include "mem.h"

typedef struct
{
    double *data;
    size_t size;
} vector_t;

vector_t *vector_alloc(const size_t n);
void vector_free(vector_t *v);

void vector_set(vector_t *v, const size_t i, double x);
double vector_get(const vector_t *v, const size_t i);

int vector_memcpy(vector_t *dst, const vector_t *src);

int vector_mul(vector_t *a, const vector_t *b);
int vector_sub(vector_t *a, const vector_t *b);
int vector_add(vector_t *a, const vector_t *b);
int vector_scale(vector_t *a, const double x);
int vector_add_constant(vector_t *a, const double x);

double vector_min(const vector_t *v);

void vector_set_all(vector_t *v, double x);
void vector_set_zero(vector_t *v);
void vector_set_ones(vector_t *v);

vector_t *vector_concat(const vector_t *v, const double x);
bool has_nonpositive_elements(const vector_t *v);

#endif  /* VECTOR_H_ */