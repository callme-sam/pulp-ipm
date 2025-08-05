#ifndef UTILS_H_
#define UTILS_H_

#include <stdbool.h>

#include "matrix.h"
#include "vector.h"

void generate_lp(matrix_t **A, vector_t **b, vector_t **c);

#endif  /* UTILS_H_ */