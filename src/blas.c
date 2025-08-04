#include <gsl/gsl_blas.h>

#include "blas.h"

int blas_dgemv(blas_trans_e trans, double alpha, const matrix_t *A, const vector_t *x, double beta, vector_t *y)
{
    if (trans == BLAS_TRANSPOSE)
        return gsl_blas_dgemv(CblasTrans, alpha, A->gsl, x->gsl, beta, y->gsl);
    else
        return gsl_blas_dgemv(CblasNoTrans, alpha, A->gsl, x->gsl, beta, y->gsl);
}

int blas_ddot(const vector_t *x, const vector_t *y, double *result)
{
    return gsl_blas_ddot(x->gsl, y->gsl, result);
}

int blas_daxpy(double alpha, const vector_t *x, vector_t *y)
{
    return gsl_blas_daxpy(alpha, x->gsl, y->gsl);
}