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

    log_matrix(LOG_LEVEL_INFO, A, "A");
    log_vector(LOG_LEVEL_INFO, b, "b");
    log_vector(LOG_LEVEL_INFO, c, "c");

    sol = solve(A, b, c);

    if (sol.status == OPTIMAL) {
        LOG_INFO("OPTIMAL SOLUTION FOUND!");
        LOG_INFO("Optimal value: %f", sol.opt_val);
        log_vector(LOG_LEVEL_INFO, sol.x_opt, "x_opt");
    }

    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(c);

    return 0;
}