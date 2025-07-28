// #include <stdio.h>

#include "utils.h"

void print_matrix(const gsl_matrix *A)
{
    size_t rows;
    size_t cols;

    rows = A->size1;
    cols = A->size2;

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            printf("%.2f ", gsl_matrix_get(A, i, j));
        }
        printf("\n");
    }
}