#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

#include "ipm.h"
#include "utils.h"

static double rand_uniform(gsl_rng *rng)
{
    return gsl_rng_uniform(rng);
}

static double rand_normal(gsl_rng *rng) {
    return gsl_ran_gaussian(rng, 1.0);
}

void generate_data(gsl_matrix **A, gsl_vector **b, gsl_vector **c)
{
    gsl_vector *tmp;
    size_t rows;
    size_t cols;

    rows = (*A)->size1;
    cols = (*A)->size2;

    assert(rows <= cols);

    tmp = gsl_vector_alloc(cols);

    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, time(NULL));

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double val = rand_normal(rng);
            gsl_matrix_set(*A, i, j, val);
        }
    }

    // A[0, :] = rand + 0.1 (boundedness)
    for (size_t j = 0; j < cols; j++) {
        double val = rand_uniform(rng) + 0.1;
        gsl_matrix_set(*A, 0, j, val);
    }

    for (size_t j = 0; j < cols; j++) {
        double val = rand_uniform(rng) + 0.01;
        gsl_vector_set(tmp, j, val);
    }

    // b = A * tmp
    gsl_blas_dgemv(CblasNoTrans, 1.0, *A, tmp, 0.0, *b);

    // c ~ U(0,1)
    for (size_t j = 0; j < cols; j++) {
        gsl_vector_set(*c, j, rand_uniform(rng));
    }

    gsl_vector_free(tmp);
    gsl_rng_free(rng);
}

int main()
{
    solution_t sol;
    size_t rows;
    size_t cols;
    
    gsl_matrix *A;
    gsl_vector *b;
    gsl_vector *c;

    rows = 5;
    cols = 5;

    A = gsl_matrix_alloc(rows, cols);
    b = gsl_vector_alloc(rows);
    c = gsl_vector_alloc(cols);

    generate_data(&A, &b, &c);

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