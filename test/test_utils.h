#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include "ipm.h"
#include "matrix.h"
#include "vector.h"

void save_problem_data(const char* filename, const matrix_t* A, const vector_t* b, const vector_t* c);
void load_problem_data(const char* filename, matrix_t** A, vector_t** b, vector_t** c);
void write_solution(const solution_t* sol, const char*sts_filename, const char* val_filename, const char* vec_filename);

#endif // TEST_UTILS_H