#include "ipm.h"

void fill_matrix(gsl_matrix *M)
{
    size_t rows;
    size_t cols;

    rows = M->size1;
    cols = M->size2;

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double val;

            val = (i * j) + j;
            gsl_matrix_set(M, i, j, val);
        }
    }
}