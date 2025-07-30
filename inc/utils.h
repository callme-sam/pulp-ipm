#ifndef UTILS_H_
#define UTILS_H_

#include <stdbool.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void fill_matrix(gsl_matrix *M);
void fill_vector(gsl_vector *v);

void print_matrix(const gsl_matrix *A);
void print_vector(const gsl_vector *v);

bool has_nonpositive_elements(const gsl_vector *v);

gsl_vector *vector_ones(size_t n);
gsl_vector *vector_concat(const gsl_vector *v, double b);

void generate_lp(gsl_matrix **A, gsl_vector **b, gsl_vector **c);

#endif  /* UTILS_H_ */