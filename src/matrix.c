#include "matrix.h"

matrix_t *matrix_alloc(const size_t n1, const size_t n2)
{
    matrix_t *m;

    m = malloc(sizeof(matrix_t));
    if (!m)
        return NULL;

    m->gsl = gsl_matrix_calloc(n1, n2);
    if (!m->gsl) {
        free(m);
        return NULL;
    }

    m->size1 = m->gsl->size1;
    m->size2 = m->gsl->size2;

    return m;
}

void matrix_free(matrix_t *m)
{
    if (!m) return;

    if (m->gsl)
        gsl_matrix_free(m->gsl);
    free(m);
}

void matrix_set(matrix_t *m, const size_t i, const size_t j, const double x)
{
    gsl_matrix_set(m->gsl, i, j, x);
}

double matrix_get(const matrix_t *m, const size_t i, const size_t j)
{
    return gsl_matrix_get(m->gsl, i, j);
}

int matrix_transpose_memcpy(matrix_t *dst, const matrix_t *src)
{
    return gsl_matrix_transpose_memcpy(dst->gsl, src->gsl);
}

int matrix_memcpy(matrix_t *dst, const matrix_t *src)
{
    return gsl_matrix_memcpy(dst->gsl, src->gsl);
}