#include <stdio.h>
#include <stdlib.h>

#include "ipm.h"
#include "utils.h"

int main()
{
    solution_t sol;
    size_t rows;
    size_t cols;
    
    gsl_matrix *A;
    gsl_vector *b;
    gsl_vector *c;

    rows = 5;
    cols = 7;

    A = gsl_matrix_alloc(rows, cols);
    b = gsl_vector_alloc(rows);
    c = gsl_vector_alloc(cols);

    generate_lp(&A, &b, &c);

    printf("A:\n");
    print_matrix(A);
    printf("b:\n");
    print_vector(b);
    printf("c:\n");
    print_vector(c);

    sol = solve(A, b, c);

    if (sol.status == OPTIMAL) {
        printf("OPTIMAL SOLUTION FOUND!\n");
        printf("Optimal value: %f\n", sol.opt_val);
        printf("Optimal solution: ");
        for (size_t i = 0; i < sol.x_opt->size; i++)
            printf("%f ", gsl_vector_get(sol.x_opt, i));
        printf("\n");

    }

    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(c);

    return 0;
}