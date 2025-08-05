#include <stdio.h>
#include <stdlib.h>

#include "test_utils.h"
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

    matrix_t* A;
    vector_t* b;
    vector_t* c;

    A = matrix_alloc(rows, cols);
    b = vector_alloc(rows);
    c = vector_alloc(cols);

    // Generate strictly feasible LP problem
    generate_lp(&A, &b, &c);

    // Save problem data
    save_problem_data(filename, A, b, c);

    // Clean up
    if (A) matrix_free(A);
    if (b) vector_free(b);
    if (c) vector_free(c);

    return EXIT_SUCCESS;
}