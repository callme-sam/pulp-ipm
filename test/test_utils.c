#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "blas.h"
#include "linalg.h"
#include "matrix.h"
#include "rng.h"
#include "vector.h"

#include "test_utils.h"

void save_problem_data(const char* filename, const matrix_t* A, const vector_t* b, const vector_t* c) {
    FILE* file = fopen(filename, "wb");
    if (!file) {
        perror("Error opening file for writing");
        return;
    }

    // Write matrix A
    size_t rows = A->size1;
    size_t cols = A->size2;
    fwrite(&rows, sizeof(size_t), 1, file);
    fwrite(&cols, sizeof(size_t), 1, file);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double val = matrix_get(A, i, j);
            fwrite(&val, sizeof(double), 1, file);
        }
    }

    // Write vector b
    size_t size_b = b->size;
    fwrite(&size_b, sizeof(size_t), 1, file);
    for (size_t i = 0; i < size_b; i++) {
        double val = vector_get(b, i);
        fwrite(&val, sizeof(double), 1, file);
    }

    // Write vector c
    size_t size_c = c->size;
    fwrite(&size_c, sizeof(size_t), 1, file);
    for (size_t i = 0; i < size_c; i++) {
        double val = vector_get(c, i);
        fwrite(&val, sizeof(double), 1, file);
    }

    fclose(file);
}

void load_problem_data(const char* filename, matrix_t** A, vector_t** b, vector_t** c) {
    FILE* file = fopen(filename, "rb");
    if (!file) {
        perror("Error opening file for reading");
        return;
    }

    // Read matrix A
    size_t rows, cols;
    int ret;
    ret = fread(&rows, sizeof(size_t), 1, file);
    if (ret != 1) return;
    ret = fread(&cols, sizeof(size_t), 1, file);
    if (ret != 1) return;
    *A = matrix_alloc(rows, cols);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double val;
            ret = fread(&val, sizeof(double), 1, file);
            if (ret != 1) return;
            matrix_set(*A, i, j, val);
        }
    }

    // Read vector b
    size_t size_b;
    ret = fread(&size_b, sizeof(size_t), 1, file);
    if (ret != 1) return;
    *b = vector_alloc(size_b);
    for (size_t i = 0; i < size_b; i++) {
        double val;
        ret = fread(&val, sizeof(double), 1, file);
        vector_set(*b, i, val);
    }

    // Read vector c
    size_t size_c;
    ret = fread(&size_c, sizeof(size_t), 1, file);
    if (ret != 1) return;
    *c = vector_alloc(size_c);
    for (size_t i = 0; i < size_c; i++) {
        double val;
        ret = fread(&val, sizeof(double), 1, file);
        if (ret != 1) return;
        vector_set(*c, i, val);
    }

    fclose(file);
}

void write_solution(const solution_t* sol, const char* sts_filename, const char* val_filename, const char* vec_filename) {

    switch (sol->status) {
        case OPTIMAL:
        {
            // Write status value
            FILE* sts_file = fopen(sts_filename, "w");
            if (sts_file) {
                fprintf(sts_file, "optimal\n");
                fclose(sts_file);
            }

            // Write optimal value
            FILE* val_file = fopen(val_filename, "w");
            if (val_file) {
                fprintf(val_file, "%.15f\n", sol->opt_val);
                fclose(val_file);
            }

            // Write solution vector
            FILE* vec_file = fopen(vec_filename, "w");
            if (vec_file && sol->x_opt) {
                for (size_t i = 0; i < sol->x_opt->size; i++) {
                    fprintf(vec_file, "%.15f\n", vector_get(sol->x_opt, i));
                }
                fclose(vec_file);
            }
            break;
        }
        case INFEASIBLE:
        {
            // Write status value
            FILE* sts_file = fopen(sts_filename, "w");
            if (sts_file) {
                fprintf(sts_file, "infeasible\n");
                fclose(sts_file);
            }
            break;
        }
        default:
        {
            // Write status value
            FILE* sts_file = fopen(sts_filename, "w");
            if (sts_file) {
                fprintf(sts_file, "failure\n");
                fclose(sts_file);
            }
        }
    }
}