#include <gsl/gsl_linalg.h>

#include "linalg.h"

int linalg_lu_decomp(matrix_t *A, permutation_t *p, int *signum)
{
    return gsl_linalg_LU_decomp(A->gsl, p->gsl, signum);
}

int linalg_lu_solve(const matrix_t *LU, const permutation_t *p, const vector_t *v, vector_t *x)
{
    return gsl_linalg_LU_solve(LU->gsl, p->gsl, v->gsl, x->gsl);
}

int linalg_sv_decomp(matrix_t *A, matrix_t *V, vector_t *S, vector_t *work)
{
    return gsl_linalg_SV_decomp(A->gsl, V->gsl, S->gsl, work->gsl);
}