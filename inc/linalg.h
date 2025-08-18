#ifndef LINALG_H_
#define LINALG_H_

#include "matrix.h"
#include "mem.h"
#include "permutation.h"
#include "vector.h"

int linalg_lu_decomp(matrix_t *A, permutation_t *p, int *signum);
int linalg_lu_solve(const matrix_t *LU, const permutation_t *p, const vector_t *v, vector_t *x);

int linalg_sv_decomp(matrix_t *A, matrix_t *V, vector_t *S);


#endif