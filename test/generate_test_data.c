#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

int main(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s rows cols output_file\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    size_t rows = atoi(argv[1]);
    size_t cols = atoi(argv[2]);
    const char* filename = argv[3];
    
    if (cols <= rows) {
        fprintf(stderr, "Error: cols must be greater than rows (cols > rows)\n");
        return EXIT_FAILURE;
    }
    
    gsl_matrix* A;
    gsl_vector* b;
    gsl_vector* c;

    A = gsl_matrix_alloc(rows, cols);
    b = gsl_vector_alloc(rows);
    c = gsl_vector_alloc(cols);
    
    // Generate strictly feasible LP problem
    generate_lp(&A, &b, &c);
    
    // Save problem data
    save_problem_data(filename, A, b, c);
    
    // Clean up
    if (A) gsl_matrix_free(A);
    if (b) gsl_vector_free(b);
    if (c) gsl_vector_free(c);
    
    return EXIT_SUCCESS;
}