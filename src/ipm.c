#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "ipm.h"
#include "utils.h"

#define MAX_ITER    (100)
#define PINV_TOL    (1e-10)
#define ALPHA       (0.25)
#define BETA        (0.5)
#define CONV        (1e-6)
#define TOL         (1e-3)
#define MU          (10.0)

const char *to_string(problem_status_e status)
{
    switch(status) {
        case OPTIMAL: return "optimal";
        case INFEASIBLE: return "infeasible";
        case FAILURE: return "failure";
        default: return "Error";
    }
}

/**
 * @brief Initialize a solution structure with default values
 *
 * @param sol Pointer to the solution structure to initialize
 */
void solution_init(solution_t *sol) {
    if (!sol) return;

    sol->status = FAILURE;
    sol->x_opt = NULL;
    sol->v_opt = NULL;
    sol->duality_gaps = NULL;
    sol->newton_steps = NULL;
    sol->opt_val = 0.0;
    sol->num_iters = 0;
}

/**
 * @brief Free all memory allocated in a solution structure
 *
 * @param sol Pointer to the solution structure to free
 */
void solution_free(solution_t *sol) {
    if (!sol) return;

    if (sol->x_opt) {
        gsl_vector_free(sol->x_opt);
        sol->x_opt = NULL;
    }

    if (sol->v_opt) {
        gsl_vector_free(sol->v_opt);
        sol->v_opt = NULL;
    }

    if (sol->duality_gaps) {
        gsl_vector_free(sol->duality_gaps);
        sol->duality_gaps = NULL;
    }

    if (sol->newton_steps) {
        gsl_vector_free(sol->newton_steps);
        sol->newton_steps = NULL;
    }

    sol->opt_val = 0.0;
    sol->num_iters = 0;
    sol->status = FAILURE;
}

/**
 * @brief Create a deep copy of a solution structure
 *
 * @param dest Destination solution (must be initialized)
 * @param src Source solution to copy from
 * @return int 0 on success, -1 on error
 */
int solution_copy(solution_t *dest, const solution_t *src) {
    if (!dest || !src) return -1;

    // Copy basic fields
    dest->status = src->status;
    dest->opt_val = src->opt_val;
    dest->num_iters = src->num_iters;

    // Copy vectors (deep copy)
    if (src->x_opt) {
        dest->x_opt = gsl_vector_alloc(src->x_opt->size);
        if (!dest->x_opt) goto error;
        gsl_vector_memcpy(dest->x_opt, src->x_opt);
    }

    if (src->v_opt) {
        dest->v_opt = gsl_vector_alloc(src->v_opt->size);
        if (!dest->v_opt) goto error;
        gsl_vector_memcpy(dest->v_opt, src->v_opt);
    }

    if (src->duality_gaps) {
        dest->duality_gaps = gsl_vector_alloc(src->duality_gaps->size);
        if (!dest->duality_gaps) goto error;
        gsl_vector_memcpy(dest->duality_gaps, src->duality_gaps);
    }

    if (src->newton_steps) {
        dest->newton_steps = gsl_vector_alloc(src->newton_steps->size);
        if (!dest->newton_steps) goto error;
        gsl_vector_memcpy(dest->newton_steps, src->newton_steps);
    }

    return 0;

error:
    solution_free(dest);
    return -1;
}

/**
 * @brief
 *
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 *
 *
 */
static void compute_residuals(const gsl_matrix *A, const gsl_vector *b, const gsl_vector *c, const gsl_vector *x,
                                const gsl_vector *v, gsl_vector *r_dual, gsl_vector *r_primal, double *r_norm)
{
    LOG_DEBUG("Computing residuals\n");

    size_t n = x->size;

    // r_dual = c - 1/x + A^T*v
    for (size_t i = 0; i < n; i++)
        gsl_vector_set(r_dual, i, gsl_vector_get(c, i) - 1.0/gsl_vector_get(x, i));
    gsl_blas_dgemv(CblasTrans, 1.0, A, v, 1.0, r_dual);

    // r_primal = A*x - b
    gsl_blas_dgemv(CblasNoTrans, 1.0, A, x, 0.0, r_primal);
    gsl_vector_sub(r_primal, b);

    // Compute ||r_dual||^2 + ||r_primal||^2
    double norm_dual, norm_primal;
    gsl_blas_ddot(r_dual, r_dual, &norm_dual);
    gsl_blas_ddot(r_primal, r_primal, &norm_primal);
    *r_norm = sqrt(norm_dual + norm_primal);
}

/**
 * @brief
 *
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 *
 *
 * @return int {description}
 */
static int compute_newton_step(const gsl_matrix *A, const gsl_vector *r_dual, const gsl_vector *r_primal,
                                const gsl_vector *x, gsl_vector *dx, gsl_vector *dv, gsl_matrix *M, gsl_vector *h_inv,
                                gsl_permutation *perm)
{
    LOG_DEBUG("Computing newton step");

    const size_t n = x->size;
    const size_t m = r_primal->size;

    // Build h_inv = x^2 (element-wise)
    for (size_t i = 0; i < n; i++) {
        double xi = gsl_vector_get(x, i);
        gsl_vector_set(h_inv, i, xi * xi);
    }

    // Compute A_hinv_rt = A * (h_inv * r_dual)
    gsl_vector *A_hinv_rt = gsl_vector_alloc(m);
    gsl_vector *temp = gsl_vector_alloc(n);
    gsl_vector_memcpy(temp, r_dual);
    gsl_vector_mul(temp, h_inv);
    gsl_blas_dgemv(CblasNoTrans, 1.0, A, temp, 0.0, A_hinv_rt);

    // Compute M = A * diag(h_inv) * A^T
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            double sum = 0.0;
            for (size_t k = 0; k < n; k++) {
                sum += gsl_matrix_get(A, i, k) *
                       gsl_vector_get(h_inv, k) *
                       gsl_matrix_get(A, j, k);
            }
            gsl_matrix_set(M, i, j, sum);
        }
    }

    // Compute rhs = r_primal - A_hinv_rt
    gsl_vector *rhs = gsl_vector_alloc(m);
    gsl_vector_memcpy(rhs, r_primal);
    gsl_vector_sub(rhs, A_hinv_rt);

    // Solve M*dv = rhs
    int signum;
    int status = gsl_linalg_LU_decomp(M, perm, &signum);
    if (status != GSL_SUCCESS) {
        LOG_ERROR("LU decomp failed: %s\n", gsl_strerror(status));
        gsl_vector_free(A_hinv_rt);
        gsl_vector_free(temp);
        gsl_vector_free(rhs);
        return -1;
    }

    status = gsl_linalg_LU_solve(M, perm, rhs, dv);
    gsl_vector_free(A_hinv_rt);
    gsl_vector_free(temp);
    gsl_vector_free(rhs);
    if (status != GSL_SUCCESS) {
        LOG_ERROR("LU solve failed: %s\n", gsl_strerror(status));
        return -1;
    }

    // Compute dx = -h_inv * (r_dual + A^T*dv)
    gsl_blas_dgemv(CblasTrans, 1.0, A, dv, 0.0, dx);  // dx = A^T*dv
    gsl_vector_add(dx, r_dual);                        // dx = r_dual + A^T*dv
    for (size_t i = 0; i < n; i++) {
        double val = -gsl_vector_get(h_inv, i) * gsl_vector_get(dx, i);
        gsl_vector_set(dx, i, val);
    }

    return 0;
}

/**
 * @brief
 *
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 *
 *
 * @return double {description}
 */
static double backtracking_line_search(const gsl_matrix *A, const gsl_vector *b, const gsl_vector *c,
                                        const gsl_vector *x, const gsl_vector *v, const gsl_vector *dx,
                                        const gsl_vector *dv, double r_norm, double alpha, double beta,
                                        gsl_vector *x_new, gsl_vector *v_new, gsl_vector *r_dual, gsl_vector *r_primal)
{
    LOG_DEBUG("Backtracking line search");
    const size_t n = x->size;
    double t = 1.0;

    // Step 1: keep x g.t. 0
    while (1) {
        bool positive = true;
        for (size_t i = 0; i < n; i++) {
            double xi = gsl_vector_get(x, i);
            double dxi = gsl_vector_get(dx, i);
            if (xi + t * dxi <= 0) {
                positive = false;
                break;
            }
        }
        if (positive) break;
        t *= beta;
    }

    // Step 2: reduce residual norm
    while (1) {
        // x_new = x + t*dx
        gsl_vector_memcpy(x_new, x);
        gsl_blas_daxpy(t, dx, x_new);

        // v_new = v + t*dv
        gsl_vector_memcpy(v_new, v);
        gsl_blas_daxpy(t, dv, v_new);

        // Compute new residuals
        double r_norm_new;
        compute_residuals(A, b, c, x_new, v_new, r_dual, r_primal, &r_norm_new);

        // Check
        if (r_norm_new <= (1 - alpha * t) * r_norm) {
            break;
        } else {
            t *= beta;
        }
    }

    return t;
}

/**
 * @brief Solves the LP centering problem using infeasible start Newton method
 *
 * @param A Constraint matrix (m x n)
 * @param b Right-hand side vector (m)
 * @param c Cost vector (n)
 * @param x0 Initial point (n)
 *
 * @return solution_e
 */
static solution_t solve_centering(const gsl_matrix *A, const gsl_vector *b, const gsl_vector *c, gsl_vector *x0)
{
    solution_t sol;
    solution_init(&sol);
    sol.status = FAILURE;

    LOG_DEBUG("Solving centering");

    const size_t m = A->size1;
    const size_t n = A->size2;
    const double convergence_threshold = 1e-6;
    const int max_iter = 100;
    const double alpha = ALPHA;
    const double beta = BETA;

    // Allocazione memoria
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector *v = gsl_vector_alloc(m);
    gsl_vector *r_dual = gsl_vector_alloc(n);
    gsl_vector *r_primal = gsl_vector_alloc(m);
    gsl_vector *dx = gsl_vector_alloc(n);
    gsl_vector *dv = gsl_vector_alloc(m);
    gsl_vector *x_new = gsl_vector_alloc(n);
    gsl_vector *v_new = gsl_vector_alloc(m);
    gsl_matrix *M = gsl_matrix_alloc(m, m);
    gsl_vector *work = gsl_vector_alloc(n);
    gsl_permutation *perm = gsl_permutation_alloc(m);

    if (!x || !v || !r_dual || !r_primal || !dx || !dv ||
        !x_new || !v_new || !M || !work || !perm) {
        goto cleanup;
    }

    // Inizializza variabili
    gsl_vector_memcpy(x, x0);
    gsl_vector_set_zero(v);

    // Loop di Newton
    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        double r_norm;

        // 1. Calcola residui
        compute_residuals(A, b, c, x, v, r_dual, r_primal, &r_norm);

        // 2. Controlla convergenza
        if (r_norm <= convergence_threshold) {
            sol.status = OPTIMAL;
            sol.x_opt = gsl_vector_alloc(n);
            sol.v_opt = gsl_vector_alloc(m);
            if (sol.x_opt && sol.v_opt) {
                gsl_vector_memcpy(sol.x_opt, x);
                gsl_vector_memcpy(sol.v_opt, v);
            }
            sol.num_iters = iter + 1;
            break;
        }

        // 3. Calcola passi di Newton
        if (compute_newton_step(A, r_dual, r_primal, x, dx, dv, M, work, perm) != 0) {
            break;
        }

        // 4. Backtracking line search
        backtracking_line_search(A, b, c, x, v, dx, dv, r_norm, alpha, beta, x_new, v_new, r_dual, r_primal);

        // 5. Aggiorna variabili
        gsl_vector_memcpy(x, x_new);
        gsl_vector_memcpy(v, v_new);
    }

    if (iter == max_iter) {
        sol.status = FAILURE;
    }

cleanup:
    // Libera memoria
    if (x) gsl_vector_free(x);
    if (v) gsl_vector_free(v);
    if (r_dual) gsl_vector_free(r_dual);
    if (r_primal) gsl_vector_free(r_primal);
    if (dx) gsl_vector_free(dx);
    if (dv) gsl_vector_free(dv);
    if (x_new) gsl_vector_free(x_new);
    if (v_new) gsl_vector_free(v_new);
    if (M) gsl_matrix_free(M);
    if (work) gsl_vector_free(work);
    if (perm) gsl_permutation_free(perm);

    return sol;
    return sol;
}

/**
 * @brief
 *
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 *
 *
 * @return solution_t {description}
 */
static solution_t solve_feasible_start(const gsl_matrix *A, const gsl_vector *b, const gsl_vector *c, const gsl_vector *x0)
{
    solution_t sol;
    gsl_vector *tc;
    gsl_vector *x;
    double gap;
    double t;

    LOG_INFO("solving feasible start");

    solution_init(&sol);
    t = 1.0;

    sol.duality_gaps = gsl_vector_alloc(MAX_ITER);
    sol.newton_steps = gsl_vector_alloc(MAX_ITER);
    tc = gsl_vector_alloc(c->size);
    x = gsl_vector_alloc(x0->size);
    gsl_vector_memcpy(x, x0);

    for (sol.num_iters = 0; sol.num_iters < MAX_ITER; sol.num_iters++) {
        solution_t center_sol;

        // Solve centering problem with current t
        gsl_vector_memcpy(tc, c);
        gsl_vector_scale(tc, t);

        center_sol = solve_centering(A, b, tc, x);
        if (center_sol.status != OPTIMAL) {
            LOG_ERROR("Centering step failed at iteration %d", sol.num_iters);
            solution_free(&center_sol);
            break;
        }

        // Update solution
        gsl_vector_memcpy(x, center_sol.x_opt);

        // Compute duality gap
        gap = A->size2 / t;
        gsl_vector_set(sol.duality_gaps, sol.num_iters, gap);
        gsl_vector_set(sol.newton_steps, sol.num_iters, center_sol.num_iters);

        // Check convergence
        if (gap < TOL) {
            sol.status = OPTIMAL;
            sol.x_opt = gsl_vector_alloc(x->size);
            sol.v_opt = gsl_vector_alloc(center_sol.v_opt->size);
            gsl_vector_memcpy(sol.x_opt, x);
            gsl_vector_memcpy(sol.v_opt, center_sol.v_opt);
            gsl_blas_ddot(c, x, &sol.opt_val);
            break;
        }

        t = t * MU;
        solution_free(&center_sol);
    }

    LOG_INFO("Found %s solution", to_string(sol.status));

    gsl_vector_free(tc);
    gsl_vector_free(x);

    return sol;
}

/**
 * @brief
 *
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 * @param {type} {name} {description}
 *
 *
 * @return problem_status_e {description}
 */
static solution_t solve_auxiliary_lp(const gsl_matrix *A, const gsl_vector *b, const gsl_vector *x)
{
    solution_t aux_sol;
    gsl_vector *A_ones;
    gsl_vector *ones;
    gsl_vector *tmp;
    solution_t sol;
    gsl_matrix *A1;
    gsl_vector *b1;
    gsl_vector *c1;
    gsl_vector *z0;
    size_t rows;
    size_t cols;
    double t;

    LOG_INFO("Solving auxiliary problem");

    solution_init(&sol);
    rows = A->size1;
    cols = A->size2;

    // Step 1: Build auxiliary problem
    A1 = gsl_matrix_alloc(rows, cols + 1);
    b1 = gsl_vector_alloc(rows);
    c1 = gsl_vector_alloc(cols + 1);
    A_ones = gsl_vector_alloc(rows);
    tmp = gsl_vector_alloc(cols);

    // Build A1 = [A, -A*ones]
    ones = vector_ones(cols);
    gsl_blas_dgemv(CblasNoTrans, -1.0, A, ones, 0.0, A_ones);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++)
            gsl_matrix_set(A1, i, j, gsl_matrix_get(A, i, j));
        gsl_matrix_set(A1, i, cols, gsl_vector_get(A_ones, i));
    }

    // Build b1 = b - A*ones
    gsl_vector_memcpy(b1, b);
    gsl_blas_daxpy(-1.0, A_ones, b1);

    // Build c1 = [0, ..., 0, 1]
    gsl_vector_set_zero(c1);
    gsl_vector_set(c1, cols, 1.0);

    // Step 2: Construct initial point

    // Compute z0 = [x + (t-1)*ones, t], with t = 2 - min_x
    t = 2.0 - gsl_vector_min(x);
    gsl_vector_memcpy(tmp, x);
    gsl_blas_daxpy(t - 1.0, ones, tmp);
    z0 = vector_concat(tmp, t);

    // Step 3: Solve LP
    aux_sol = solve_feasible_start(A1, b1, c1, z0);

    // Step 4: Process results
    if (aux_sol.status == OPTIMAL && aux_sol.opt_val < 1.0) {
        LOG_INFO("Auxiliary Problem il feasible, extracting solution");
        sol.status = OPTIMAL;

        // Extract feasible point: x = z_opt[:n] - (z_opt[n] - 1)*ones
        sol.x_opt = gsl_vector_alloc(cols);
        double t_opt = gsl_vector_get(aux_sol.x_opt, cols);
        for (size_t i = 0; i < cols; i++) {
            double v = gsl_vector_get(aux_sol.x_opt, i) - (t_opt - 1.0);
            gsl_vector_set(sol.x_opt, i, v);
        }
    } else {
        LOG_WARNING("Auxiliary Problem is INFEASIBLE - status=%d - opt_val=%f - n_iters=%d", sol.status, sol.opt_val, sol.num_iters);
        sol.status = INFEASIBLE;
    }

    if (A_ones) gsl_vector_free(A_ones);
    if (ones) gsl_vector_free(ones);
    if (tmp) gsl_vector_free(tmp);
    if (A1) gsl_matrix_free(A1);
    if (b1) gsl_vector_free(b1);
    if (c1) gsl_vector_free(c1);
    if (z0) gsl_vector_free(z0);
    solution_free(&aux_sol);

    return sol;
}

/**
 * @brief Solves a linear system using the Moore–Penrose pseudo-inverse via SVD.
 *
 * Computes the minimum-norm least-squares solution to the linear system \( Ax = b \)
 * using the pseudo-inverse \( A^+ \) obtained via Singular Value Decomposition (SVD).
 * This method works for both overdetermined and underdetermined systems.
 *
 * @param[in] A Pointer to a GSL matrix of size M×N representing the system matrix.
 * @param[in] b Pointer to a GSL vector of size M representing the right-hand side.
 *
 * @return On success: a newly allocated GSL vector of size N containing the solution \( x = A^+ b \).
 *         On failure: NULL (e.g., due to memory allocation failure or invalid input).
 *
 * @warning The caller is responsible for freeing the returned vector using `gsl_vector_free()`.
 * @note The pseudo-inverse is computed via SVD as \( A^+ = V \Sigma^+ U^T \) or \( A^+ = U \Sigma^+ V^T \)
 *       depending on the shape of A.
 */
static gsl_vector *pinv_solve(const gsl_matrix *A, const gsl_vector *b) {
    LOG_INFO("Solving pseudo-inverse via SVD");

    size_t rows;
    size_t cols;

    gsl_matrix *U;
    gsl_matrix *V;

    gsl_vector *work;
    gsl_vector *tmp;
    gsl_vector *x;
    gsl_vector *S;

    rows = A->size1;
    cols = A->size2;

    if (rows >= cols) {
        U = gsl_matrix_alloc(rows, cols);
        V = gsl_matrix_alloc(cols, cols);
        work = gsl_vector_alloc(cols);
        tmp = gsl_vector_alloc(cols);
        S = gsl_vector_alloc(cols);
        x = gsl_vector_alloc(cols);

        gsl_matrix_memcpy(U, A);
        gsl_linalg_SV_decomp(U, V, S, work);

        gsl_blas_dgemv(CblasTrans, 1.0, U, b, 0.0, tmp);    // tmp = Uᵗ b

        // Pseudoinverse: tmp = S⁺ * tmp
        for (size_t i = 0; i < cols; i++) {
            double s_val = gsl_vector_get(S, i);
            if (s_val > PINV_TOL) {
                gsl_vector_set(tmp, i, gsl_vector_get(tmp, i) / s_val);
            } else {
                gsl_vector_set(tmp, i, 0.0);
            }
        }

        gsl_blas_dgemv(CblasNoTrans, 1.0, V, tmp, 0.0, x);  // x = V * Σ⁺ * tmp
    } else {
        // Transpose A and compute SVD of Aᵗ
        U = gsl_matrix_alloc(cols, rows);
        V = gsl_matrix_alloc(rows, rows);
        work = gsl_vector_alloc(rows);
        tmp = gsl_vector_alloc(rows);
        S = gsl_vector_alloc(rows);
        x = gsl_vector_alloc(cols);

        gsl_matrix_transpose_memcpy(U, A);                  // U = Aᵗ
        gsl_linalg_SV_decomp(U, V, S, work);

        gsl_blas_dgemv(CblasTrans, 1.0, V, b, 0.0, tmp);    // tmp = Vᵗ b

        // Pseudoinverse: tmp = S⁺ * tmp
        for (size_t i = 0; i < rows; i++) {
            double s_val = gsl_vector_get(S, i);
            if (s_val > PINV_TOL) {
                gsl_vector_set(tmp, i, gsl_vector_get(tmp, i) / s_val);
            } else {
                gsl_vector_set(tmp, i, 0.0);
            }
        }

        gsl_blas_dgemv(CblasNoTrans, 1.0, U, tmp, 0.0, x);  // x = V * Σ⁺ * tmp
    }

    if (work) gsl_vector_free(work);
    if (tmp) gsl_vector_free(tmp);
    if (U) gsl_matrix_free(U);
    if (V) gsl_matrix_free(V);
    if (S) gsl_vector_free(S);

    return x;
}

/**
 * @brief Performs Phase I of the interior point method to find a feasible starting point.
 *
 * This function constructs and solves an auxiliary (Phase I) linear program to determine whether the original LP
 * problem has at least one feasible point (i.e., a point x such that Ax = b and x > 0). If the problem is strictly
 * feasible, the function attempts to use the pseudo-inverse to find a starting point directly. Otherwise, it formulates
 * and solves an auxiliary LP to find one.
 *
 * @param [in] A Constraint matrix of size m x n (rows × columns).
 * @param [in] b Right-hand side vector of size m.
 *
 * @return solution_t A solution object containing:
 *         - status: OPTIMAL / INFEASIBLE / FAILURE
 *         - x_opt: Feasible starting point (if status == OPTIMAL)
 *         - Diagnostics fields may be empty
 */
static solution_t phase_one(const gsl_matrix *A, const gsl_vector *b)
{
    LOG_INFO("Solving phase I");

    solution_t sol;
    solution_init(&sol);

    // Step 1: Compute initial point via pseudo-inverse
    gsl_vector *x;
    x = pinv_solve(A, b);
    if (!x) {
        LOG_ERROR("Computation of feasible starting point FAILED");
        sol.status = FAILURE;
        return sol;
    }

    // Step 2: Check feasibility
    if (has_nonpositive_elements(x)) {
        LOG_INFO("X is NOT a strictly feasible starting point, constructing Auxiliary LP Problem");

        // Step 3: Solve auxiliary LP
        solution_t aux_sol;
        aux_sol = solve_auxiliary_lp(A, b, x);
        solution_copy(&sol, &aux_sol);
        solution_free(&aux_sol);
    } else {
        LOG_INFO("X is a strictly feasible starting point");

        sol.status = OPTIMAL;
        sol.x_opt = gsl_vector_alloc(x->size);
        gsl_vector_memcpy(sol.x_opt, x);
    }

    if (x) gsl_vector_free(x);
    return sol;
}

/**
 * @brief Solves a linear programming (LP) problem using a two-phase method.
 *
 * This function applies a two-phase approach:
 * - **Phase I:** finds a feasible starting point.
 * - **Phase II:** performs optimization starting from the feasible solution found.
 *
 * If Phase I fails, the problem is declared INFEASIBLE.
 *
 * @param A Constraint matrix (m x n) of the LP problem.
 * @param b Right-hand side vector (m x 1) of the constraints.
 * @param c Coefficient vector of the objective function (n x 1).
 *
 * @return solution_t A `solution_t` object containing:
 *  - `status`: the solution status (e.g., OPTIMAL, INFEASIBLE).
 *  - `x_opt`: the optimal solution found (if any).
 *  - `value`: the optimal objective function value (if computed).
 */
solution_t solve(gsl_matrix *A, gsl_vector *b, gsl_vector *c)
{
    solution_t phase1_sol;
    solution_t phase2_sol;
    solution_t sol;

    LOG_INFO("Solving LP problem");

    phase1_sol = phase_one(A, b);
    if (phase1_sol.status == OPTIMAL) {
        LOG_INFO("Phase I got a feasible starting point, computing Phase II");
        phase2_sol = solve_feasible_start(A, b, c, phase1_sol.x_opt);
        solution_copy(&sol, &phase2_sol);
    } else {
        LOG_ERROR("Phase I failed. Original LP is INFEASIBLE");
        solution_copy(&sol, &phase1_sol);
    }

    solution_free(&phase1_sol);
    solution_free(&phase2_sol);
    return sol;
}