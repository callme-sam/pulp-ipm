#include <stdio.h>
#include <stdlib.h>

#include "ipm.h"
#include "log.h"
#include "utils.h"

#include "pmsis.h"

static int test_entry()
{
    solution_t sol;
    size_t rows;
    size_t cols;

    matrix_t *A;
    vector_t *b;
    vector_t *c;

    rows = 5;
    cols = 7;

    A = matrix_alloc(rows, cols);
    b = vector_alloc(rows);
    c = vector_alloc(cols);

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

    matrix_free(A);
    vector_free(b);
    vector_free(c);

    return 0;
}

static void test_kickoff(void *arg)
{
    int ret;

    ret = test_entry();
    pmsis_exit(ret);
}

int main()
{
    return pmsis_kickoff((void *) test_kickoff);
}