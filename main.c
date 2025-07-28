#include <stdio.h>
#include <stdlib.h>

#include "ipm.h"
#include "utils.h"

int main()
{
    int m = 2;
    int n = 3;

    gsl_matrix *A = gsl_matrix_alloc(m, n);

    fill_matrix(A);
    print_matrix(A);

    return 0;
}