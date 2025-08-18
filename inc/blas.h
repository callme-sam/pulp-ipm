#ifndef BLAS_H_
#define BLAS_H_

#include "matrix.h"
#include "mem.h"
#include "vector.h"

typedef enum
{
    BLAS_TRANSPOSE,
    BLAS_NO_TRANSPOSE
} blas_trans_e;

int blas_dgemv(blas_trans_e trans, double alpha, const matrix_t *A, const vector_t *x, double beta, vector_t *y);
int blas_ddot(const vector_t *x, const vector_t *y, double *result);
int blas_daxpy(double alpha, const vector_t *x, vector_t *y);

#endif  /* BLAS_H_ */