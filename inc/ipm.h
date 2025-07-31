#ifndef IPM_H_
#define IPM_H_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef enum {
    OPTIMAL     = 0,
    INFEASIBLE  = 1,
    FAILURE     = 2,
    NOT_SET     = 3,
} problem_status_e;

typedef struct {
    problem_status_e status;
    gsl_vector *duality_gaps;   // History of duality gaps
    gsl_vector *newton_steps;   // History of newton steps
    gsl_vector *x_opt;          // Primal Optimal Solution
    gsl_vector *v_opt;          // Dual Optimal Solution
    double opt_val;             // Optimal Value
    int num_iters;              // Total Number of iterations
} solution_t;

const char *to_string(problem_status_e status);

void solution_init(solution_t *sol);
void solution_free(solution_t *sol);
int solution_copy(solution_t *dest, const solution_t *src);

solution_t solve(gsl_matrix *A, gsl_vector *b, gsl_vector *c);

#endif  /* IPM_H_ */