#ifndef IPM_H_
#define IPM_H_

#include "matrix.h"
#include "vector.h"

typedef enum {
    OPTIMAL     = 0,
    INFEASIBLE  = 1,
    FAILURE     = 2,
    NOT_SET     = 3,
} problem_status_e;

typedef struct {
    problem_status_e status;
    vector_t *duality_gaps;   // History of duality gaps
    vector_t *newton_steps;   // History of newton steps
    vector_t *x_opt;          // Primal Optimal Solution
    vector_t *v_opt;          // Dual Optimal Solution
    double opt_val;             // Optimal Value
    int num_iters;              // Total Number of iterations
} solution_t;

const char *to_string(problem_status_e status);

void solution_init(solution_t *sol);
void solution_free(solution_t *sol);
int solution_copy(solution_t *dest, const solution_t *src);

solution_t solve(matrix_t *A, vector_t *b, vector_t *c);

#endif  /* IPM_H_ */