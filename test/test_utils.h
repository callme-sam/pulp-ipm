#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "../inc/ipm.h"

void save_problem_data(const char* filename, const gsl_matrix* A, const gsl_vector* b, const gsl_vector* c);
void load_problem_data(const char* filename, gsl_matrix** A, gsl_vector** b, gsl_vector** c);
void write_solution(const solution_t* sol, const char*sts_filename, const char* val_filename, const char* vec_filename);

#endif // TEST_UTILS_H