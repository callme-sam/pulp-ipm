#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s m n output_file\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    size_t m = atoi(argv[1]);
    size_t n = atoi(argv[2]);
    const char* filename = argv[3];
    
    // if (n <= m) {
    //     fprintf(stderr, "Error: n must be greater than m (n > m)\n");
    //     return EXIT_FAILURE;
    // }
    
    gsl_matrix* A = NULL;
    gsl_vector* b = NULL;
    gsl_vector* c = NULL;
    
    // Generate strictly feasible LP problem
    generate_strictly_feasible_lp(&A, &b, &c, m, n);
    
    // Save problem data
    save_problem_data(filename, A, b, c);
    
    // Clean up
    if (A) gsl_matrix_free(A);
    if (b) gsl_vector_free(b);
    if (c) gsl_vector_free(c);
    
    return EXIT_SUCCESS;
}