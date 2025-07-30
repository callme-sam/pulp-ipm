#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>

#include "test_utils.h"

void save_problem_data(const char* filename, const gsl_matrix* A, const gsl_vector* b, const gsl_vector* c) {
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
            double val = gsl_matrix_get(A, i, j);
            fwrite(&val, sizeof(double), 1, file);
        }
    }
    
    // Write vector b
    size_t size_b = b->size;
    fwrite(&size_b, sizeof(size_t), 1, file);
    for (size_t i = 0; i < size_b; i++) {
        double val = gsl_vector_get(b, i);
        fwrite(&val, sizeof(double), 1, file);
    }
    
    // Write vector c
    size_t size_c = c->size;
    fwrite(&size_c, sizeof(size_t), 1, file);
    for (size_t i = 0; i < size_c; i++) {
        double val = gsl_vector_get(c, i);
        fwrite(&val, sizeof(double), 1, file);
    }
    
    fclose(file);
}

void load_problem_data(const char* filename, gsl_matrix** A, gsl_vector** b, gsl_vector** c) {
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
    *A = gsl_matrix_alloc(rows, cols);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double val;
            ret = fread(&val, sizeof(double), 1, file);
            if (ret != 1) return;
            gsl_matrix_set(*A, i, j, val);
        }
    }
    
    // Read vector b
    size_t size_b;
    ret = fread(&size_b, sizeof(size_t), 1, file);
    if (ret != 1) return;
    *b = gsl_vector_alloc(size_b);
    for (size_t i = 0; i < size_b; i++) {
        double val;
        ret = fread(&val, sizeof(double), 1, file);
        gsl_vector_set(*b, i, val);
    }
    
    // Read vector c
    size_t size_c;
    ret = fread(&size_c, sizeof(size_t), 1, file);
    if (ret != 1) return;
    *c = gsl_vector_alloc(size_c);
    for (size_t i = 0; i < size_c; i++) {
        double val;
        ret = fread(&val, sizeof(double), 1, file);
        if (ret != 1) return;
        gsl_vector_set(*c, i, val);
    }
    
    fclose(file);
}

void write_solution(const solution_t* sol, const char* val_filename, const char* vec_filename) {
    if (sol->status == OPTIMAL) {
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
                fprintf(vec_file, "%.15f\n", gsl_vector_get(sol->x_opt, i));
            }
            fclose(vec_file);
        }
    }
}