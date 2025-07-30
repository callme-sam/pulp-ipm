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

int get_matrix_rank(const gsl_matrix *A, double tol)
{
    size_t m = A->size1;
    size_t n = A->size2;

    gsl_matrix *A_work = NULL;

    if (m < n) {
        // Trasponi A se ha più colonne che righe
        A_work = gsl_matrix_alloc(n, m);
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                gsl_matrix_set(A_work, j, i, gsl_matrix_get(A, i, j));
        size_t tmp = m;
        m = n;
        n = tmp;
    } else {
        A_work = gsl_matrix_alloc(m, n);
        gsl_matrix_memcpy(A_work, A);
    }

    gsl_vector *S = gsl_vector_alloc(n);
    gsl_vector *work = gsl_vector_alloc(n);
    gsl_matrix *V = gsl_matrix_alloc(n, n);

    int status = gsl_linalg_SV_decomp(A_work, V, S, work);
    if (status != GSL_SUCCESS) {
        fprintf(stderr, "SVD failed: %s\n", gsl_strerror(status));
        gsl_matrix_free(A_work);
        gsl_vector_free(S);
        gsl_vector_free(work);
        gsl_matrix_free(V);
        return -1;
    }

    double max_singular = gsl_vector_get(S, 0);
    for (size_t i = 1; i < n; ++i) {
        double s = gsl_vector_get(S, i);
        if (s > max_singular)
            max_singular = s;
    }

    if (tol < 0)
        tol = max_singular * GSL_DBL_EPSILON * (double)(m > n ? m : n);

    int rank = 0;
    for (size_t i = 0; i < n; ++i) {
        if (gsl_vector_get(S, i) > tol)
            rank++;
    }

    gsl_matrix_free(A_work);
    gsl_vector_free(S);
    gsl_vector_free(work);
    gsl_matrix_free(V);

    return rank;
}

void generate_strictly_feasible_lp(gsl_matrix** A, gsl_vector** b, gsl_vector** c, size_t m, size_t n) {
    // Initialize random number generator
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);
    
    // Allocate memory for the problem
    *A = gsl_matrix_alloc(m, n);
    *b = gsl_vector_alloc(m);
    *c = gsl_vector_alloc(n);
    
    // Generate random strictly feasible solution x0 > 0
    gsl_vector* x0 = gsl_vector_alloc(n);
    for (size_t i = 0; i < n; i++) {
        gsl_vector_set(x0, i, gsl_rng_uniform(r) + 0.1); // Ensure x0 > 0
    }
    
    // Generate random cost vector c
    for (size_t i = 0; i < n; i++) {
        gsl_vector_set(*c, i, gsl_ran_gaussian(r, 1.0));
    }
    
    // Generate random matrix A with full row rank
    do {
        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < n; j++) {
                gsl_matrix_set(*A, i, j, gsl_ran_gaussian(r, 1.0));
            }
        }
    } while (get_matrix_rank(*A, -1.0) < m); // Ensure full row rank
    
    // Compute b = A * x0 to ensure feasibility
    gsl_blas_dgemv(CblasNoTrans, 1.0, *A, x0, 0.0, *b);
    
    // Cleanup
    gsl_vector_free(x0);
    gsl_rng_free(r);
}

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
    (void) fread(&rows, sizeof(size_t), 1, file);
    (void) fread(&cols, sizeof(size_t), 1, file);
    *A = gsl_matrix_alloc(rows, cols);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double val;
            (void) fread(&val, sizeof(double), 1, file);
            gsl_matrix_set(*A, i, j, val);
        }
    }
    
    // Read vector b
    size_t size_b;
    (void) fread(&size_b, sizeof(size_t), 1, file);
    *b = gsl_vector_alloc(size_b);
    for (size_t i = 0; i < size_b; i++) {
        double val;
        (void) fread(&val, sizeof(double), 1, file);
        gsl_vector_set(*b, i, val);
    }
    
    // Read vector c
    size_t size_c;
    (void) fread(&size_c, sizeof(size_t), 1, file);
    *c = gsl_vector_alloc(size_c);
    for (size_t i = 0; i < size_c; i++) {
        double val;
        (void) fread(&val, sizeof(double), 1, file);
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