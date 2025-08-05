#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdlib.h>

#include <gsl/gsl_matrix.h>

typedef struct
{
    double *data;
    size_t size1;
    size_t size2;

    // TODO: remove
    gsl_matrix *gsl;
} matrix_t;

matrix_t *matrix_alloc(const size_t n1, const size_t n2);
void matrix_free(matrix_t *m);

void matrix_set(matrix_t *m, const size_t i, const size_t j, const double x);
double matrix_get(const matrix_t *m, const size_t i, const size_t j);

int matrix_transpose_memcpy(matrix_t *dst, const matrix_t *src);
int matrix_memcpy(matrix_t *dst, const matrix_t *src);

#endif  /* MATRIX_H_ */