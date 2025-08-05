#include <stdio.h>
#include <stdlib.h>

#include "ipm.h"
#include "test_utils.h"

int main(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <input_file> <output_status_file> <output_val_file> <output_vec_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    matrix_t* A = NULL;
    vector_t* b = NULL;
    vector_t* c = NULL;

    // Load problem data
    load_problem_data(argv[1], &A, &b, &c);

    if (!A || !b || !c) {
        fprintf(stderr, "Error loading problem data\n");
        return EXIT_FAILURE;
    }

    // Solve the LP problem
    solution_t sol = solve(A, b, c);

    // Write solution
    write_solution(&sol, argv[2], argv[3], argv[4]);

    // Clean up
    solution_free(&sol);
    if (A) matrix_free(A);
    if (b) vector_free(b);
    if (c) vector_free(c);

    return EXIT_SUCCESS;
}