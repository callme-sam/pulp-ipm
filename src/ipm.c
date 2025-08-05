#include <math.h>

#include "blas.h"
#include "linalg.h"
#include "log.h"

#include "ipm.h"

#define MAX_ITER    (100)
#define PINV_TOL    (1e-10)
#define ALPHA       (0.25)
#define BETA        (0.5)
#define CONV        (1e-6)
#define TOL         (1e-3)
#define MU          (10.0)

/**
 * @brief Converts a problem status code to its string representation.
 *
 * This function maps an enumerated value of type `problem_status_e` to a
 * human-readable string. Useful for logging or displaying the status of a
 * linear programming solution.
 *
 * @param status The problem status to convert (e.g., OPTIMAL, INFEASIBLE, FAILURE).
 *
 * @return A string corresponding to the given status:
 *         - "optimal" if status == OPTIMAL
 *         - "infeasible" if status == INFEASIBLE
 *         - "failure" if status == FAILURE
 *         - "Error" for any unrecognized status
 */
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
 * @brief Initializes a solution structure with default values.
 *
 * This function sets the fields of a `solution_t` structure to safe default values:
 * - `status` is set to `FAILURE`
 * - `x_opt`, `v_opt`, `duality_gaps`, and `newton_steps` are set to NULL
 * - `opt_val` is set to 0.0
 * - `num_iters` is set to 0
 *
 * This should be called before the structure is used, to avoid undefined behavior.
 *
 * @param[in,out] sol Pointer to the `solution_t` structure to initialize.
 *                    If NULL, the function does nothing.
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
 * @brief Frees all dynamically allocated memory in a solution structure.
 *
 * This function deallocates any memory previously allocated for the fields of
 * a `solution_t` structure, including:
 * - `x_opt`: optimal primal solution vector
 * - `v_opt`: optimal dual variables vector
 * - `duality_gaps`: vector of duality gap values at each iteration
 * - `newton_steps`: vector of Newton step counts at each iteration
 *
 * After deallocation, all pointers are set to `NULL`, and numeric fields are reset
 * to their default values (`opt_val = 0.0`, `num_iters = 0`, `status = FAILURE`).
 *
 * It is safe to call this function multiple times on the same structure.
 *
 * @param[in,out] sol Pointer to the `solution_t` structure to clean up.
 *                    If NULL, the function does nothing.
 */
void solution_free(solution_t *sol) {
    if (!sol) return;

    if (sol->x_opt) {
        vector_free(sol->x_opt);
        sol->x_opt = NULL;
    }

    if (sol->v_opt) {
        vector_free(sol->v_opt);
        sol->v_opt = NULL;
    }

    if (sol->duality_gaps) {
        vector_free(sol->duality_gaps);
        sol->duality_gaps = NULL;
    }

    if (sol->newton_steps) {
        vector_free(sol->newton_steps);
        sol->newton_steps = NULL;
    }

    sol->opt_val = 0.0;
    sol->num_iters = 0;
    sol->status = FAILURE;
}

/**
 * @brief Creates a deep copy of a solution structure.
 *
 * This function copies the contents of a `solution_t` structure from a source to a destination,
 * including all dynamically allocated vectors. The destination structure (`dest`) must be
 * properly initialized (e.g., via `solution_init`) before calling this function, and must not
 * contain already allocated memory to avoid leaks.
 *
 * If any memory allocation fails during the copy process, the destination structure is cleaned
 * up with `solution_free()` and the function returns -1.
 *
 * @param[out] dest Pointer to the destination `solution_t` structure.
 * @param[in]  src  Pointer to the source `solution_t` structure to copy.
 *
 * @return 0 on success, -1 on failure (e.g. memory allocation error or NULL pointers).
 *
 * @note The destination is overwritten; ensure it is not holding previously allocated data.
 *       Call `solution_free(dest)` first if needed.
 */
int solution_copy(solution_t *dest, const solution_t *src) {
    if (!dest || !src) return -1;

    // Copy basic fields
    dest->status = src->status;
    dest->opt_val = src->opt_val;
    dest->num_iters = src->num_iters;

    // Copy vectors (deep copy)
    if (src->x_opt) {
        dest->x_opt = vector_alloc(src->x_opt->size);
        if (!dest->x_opt) goto error;
        vector_memcpy(dest->x_opt, src->x_opt);
    }

    if (src->v_opt) {
        dest->v_opt = vector_alloc(src->v_opt->size);
        if (!dest->v_opt) goto error;
        vector_memcpy(dest->v_opt, src->v_opt);
    }

    if (src->duality_gaps) {
        dest->duality_gaps = vector_alloc(src->duality_gaps->size);
        if (!dest->duality_gaps) goto error;
        vector_memcpy(dest->duality_gaps, src->duality_gaps);
    }

    if (src->newton_steps) {
        dest->newton_steps = vector_alloc(src->newton_steps->size);
        if (!dest->newton_steps) goto error;
        vector_memcpy(dest->newton_steps, src->newton_steps);
    }

    return 0;

error:
    solution_free(dest);
    return -1;
}

/**
 * @brief Computes the residuals and their combined norm for the Newton step.
 *
 * This function computes the dual and primal residuals used in the infeasible start Newton method
 * for solving linear programming problems in barrier form:
 *
 * - Dual residual: \f$ r_{\text{dual}} = c - X^{-1}\mathbf{1} + A^T v \f$
 * - Primal residual: \f$ r_{\text{primal}} = Ax - b \f$
 * - Combined residual norm: \f$ \|r_{\text{dual}}\|^2 + \|r_{\text{primal}}\|^2 \f$
 *
 * These quantities are used to determine convergence and drive the Newton direction updates.
 *
 * @param[in]  A        Constraint matrix of size (m x n).
 * @param[in]  b        Right-hand side vector of size (m).
 * @param[in]  c        Cost vector of size (n).
 * @param[in]  x        Current primal variable vector (must be strictly positive).
 * @param[in]  v        Current dual variable vector.
 * @param[out] r_dual   Output dual residual vector (size n).
 * @param[out] r_primal Output primal residual vector (size m).
 * @param[out] r_norm   Output scalar for the combined residual norm.
 */
static void compute_residuals(const matrix_t *A, const vector_t *b, const vector_t *c, const vector_t *x,
                                const vector_t *v, vector_t *r_dual, vector_t *r_primal, double *r_norm)
{
    LOG_DEBUG("Computing residuals\n");

    size_t n = x->size;

    // r_dual = c - 1/x + A^T*v
    for (size_t i = 0; i < n; i++)
        vector_set(r_dual, i, vector_get(c, i) - 1.0/vector_get(x, i));
    blas_dgemv(BLAS_TRANSPOSE, 1.0, A, v, 1.0, r_dual);

    // r_primal = A*x - b
    blas_dgemv(BLAS_NO_TRANSPOSE, 1.0, A, x, 0.0, r_primal);
    vector_sub(r_primal, b);

    // Compute ||r_dual||^2 + ||r_primal||^2
    double norm_dual, norm_primal;
    blas_ddot(r_dual, r_dual, &norm_dual);
    blas_ddot(r_primal, r_primal, &norm_primal);
    *r_norm = sqrt(norm_dual + norm_primal);
}

/**
 * @brief Compute element-wise square of vector x and store in out.
 *
 * @param[in]  x   Input vector of size n.
 * @param[out] out Output vector of size n where squared values are stored.
 */
static void vector_square(const vector_t *x, vector_t *out) {
    size_t v_len;

    v_len = x->size;
    for (size_t i = 0; i < v_len; i++) {
        double xi = vector_get(x, i);
        vector_set(out, i, xi * xi);
    }
}

/**
 * @brief Compute matrix \( M = A \cdot \mathrm{diag}(h_{\text{inv}}) \cdot A^T \).
 *
 * @param[in]  A      Constraint matrix of size m×n.
 * @param[in]  h_inv  Vector of length n representing the diagonal elements.
 * @param[out] M      Output matrix of size m×m.
 */
static void compute_M_matrix(const matrix_t *A, const vector_t *h_inv, matrix_t *M) {
    size_t rows;
    size_t cols;

    rows = A->size1;
    cols = A->size2;

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < rows; j++) {
            double sum = 0.0;
            for (size_t k = 0; k < cols; k++) {
                sum += matrix_get(A, i, k) * vector_get(h_inv, k) * matrix_get(A, j, k);
            }
            matrix_set(M, i, j, sum);
        }
    }
}

/**
 * @brief Compute the Newton step for the centering problem.
 *
 * @param[in]  A       Constraint matrix (m×n).
 * @param[in]  r_dual  Dual residual vector (size n).
 * @param[in]  r_primal Primal residual vector (size m).
 * @param[in]  x       Current iterate vector (size n).
 * @param[out] dx      Newton step for primal variables (size n).
 * @param[out] dv      Newton step for dual variables (size m).
 * @param[in,out] M    Workspace matrix (m×m).
 * @param[in,out] h_inv Workspace vector (size n), stores element-wise \( x^2 \).
 * @param[in,out] perm Workspace permutation for LU decomposition (size m).
 *
 * @return int 0 on success, -1 on failure.
 */
static int compute_newton_step(const matrix_t *A, const vector_t *r_dual, const vector_t *r_primal,
                              const vector_t *x, vector_t *dx, vector_t *dv, matrix_t *M, vector_t *h_inv,
                              permutation_t *perm) {
    LOG_DEBUG("Computing Newton step");

    size_t n = x->size;
    size_t m = r_primal->size;

    vector_square(x, h_inv);

    vector_t *temp = vector_alloc(n);
    vector_t *A_hinv_rt = vector_alloc(m);
    vector_t *rhs = vector_alloc(m);

    // temp = r_dual .* h_inv (element-wise)
    vector_memcpy(temp, r_dual);
    vector_mul(temp, h_inv);

    // A_hinv_rt = A * temp
    blas_dgemv(BLAS_NO_TRANSPOSE, 1.0, A, temp, 0.0, A_hinv_rt);

    compute_M_matrix(A, h_inv, M);

    // rhs = r_primal - A_hinv_rt
    vector_memcpy(rhs, r_primal);
    vector_sub(rhs, A_hinv_rt);

    // Solve M*dv = rhs
    int signum, status = linalg_lu_decomp(M, perm, &signum);
    if (status != 0) {
        LOG_ERROR("LU decomposition failed: %s", err_to_str(status));
        goto cleanup_error;
    }
    status = linalg_lu_solve(M, perm, rhs, dv);
    if (status != 0) {
        LOG_ERROR("LU solve failed: %s", err_to_str(status));
        goto cleanup_error;
    }

    // dx = -h_inv * (r_dual + A^T * dv)
    blas_dgemv(BLAS_TRANSPOSE, 1.0, A, dv, 0.0, dx); // dx = A^T * dv
    vector_add(dx, r_dual);                         // dx += r_dual
    for (size_t i = 0; i < n; i++) {
        double val = -vector_get(h_inv, i) * vector_get(dx, i);
        vector_set(dx, i, val);
    }

    vector_free(temp);
    vector_free(A_hinv_rt);
    vector_free(rhs);
    return 0;

cleanup_error:
    vector_free(temp);
    vector_free(A_hinv_rt);
    vector_free(rhs);
    return -1;
}


/**
 * @brief Perform backtracking line search to find step size \(t\) that maintains positivity and reduces residual norm.
 *
 * @param[in]  A        Constraint matrix (m×n).
 * @param[in]  b        Right-hand side vector (size m).
 * @param[in]  c        Cost vector (size n).
 * @param[in]  x        Current primal vector (size n).
 * @param[in]  v        Current dual vector (size m).
 * @param[in]  dx       Newton step for primal variables (size n).
 * @param[in]  dv       Newton step for dual variables (size m).
 * @param[in]  r_norm   Current residual norm.
 * @param[out] x_new    Output vector for updated primal variables (size n).
 * @param[out] v_new    Output vector for updated dual variables (size m).
 * @param[in,out] r_dual Updated dual residual vector (size n).
 * @param[in,out] r_primal Updated primal residual vector (size m).
 *
 * @return double The step size \(t\) found by backtracking.
 */
static double backtracking_line_search(const matrix_t *A, const vector_t *b, const vector_t *c,
                                        const vector_t *x, const vector_t *v, const vector_t *dx,
                                        const vector_t *dv, double r_norm, vector_t *x_new, vector_t *v_new,
                                        vector_t *r_dual, vector_t *r_primal)
{
    LOG_DEBUG("Backtracking line search");
    const size_t n = x->size;
    double t = 1.0;

    // Step 1: keep x >= 0
    while (1) {
        bool positive = true;
        for (size_t i = 0; i < n; i++) {
            double xi = vector_get(x, i);
            double dxi = vector_get(dx, i);
            if (xi + t * dxi <= 0) {
                positive = false;
                break;
            }
        }
        if (positive) break;
        t *= BETA;
    }

    // Step 2: reduce residual norm
    while (1) {
        // x_new = x + t*dx
        vector_memcpy(x_new, x);
        blas_daxpy(t, dx, x_new);

        // v_new = v + t*dv
        vector_memcpy(v_new, v);
        blas_daxpy(t, dv, v_new);

        // Compute new residuals
        double r_norm_new;
        compute_residuals(A, b, c, x_new, v_new, r_dual, r_primal, &r_norm_new);

        // Check
        if (r_norm_new <= (1 - ALPHA * t) * r_norm) {
            break;
        } else {
            t *= BETA;
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
static solution_t solve_centering(const matrix_t *A, const vector_t *b, const vector_t *c, vector_t *x0)
{
    LOG_DEBUG("Solving centering");

    solution_t sol;
    solution_init(&sol);

    const size_t rows = A->size1;
    const size_t cols = A->size2;
    int iter;

    vector_t *x = vector_alloc(cols);
    vector_t *v = vector_alloc(rows);
    vector_t *r_dual = vector_alloc(cols);
    vector_t *r_primal = vector_alloc(rows);
    vector_t *dx = vector_alloc(cols);
    vector_t *dv = vector_alloc(rows);
    vector_t *x_new = vector_alloc(cols);
    vector_t *v_new = vector_alloc(rows);
    matrix_t *M = matrix_alloc(rows, rows);
    vector_t *work = vector_alloc(cols);
    permutation_t *perm = permutation_alloc(rows);

    vector_memcpy(x, x0);
    vector_set_zero(v);

    for (iter = 0; iter < MAX_ITER; iter++) {
        double r_norm;

        compute_residuals(A, b, c, x, v, r_dual, r_primal, &r_norm);

        if (r_norm <= CONV) {
            sol.status = OPTIMAL;
            sol.x_opt = vector_alloc(cols);
            sol.v_opt = vector_alloc(rows);
            if (sol.x_opt && sol.v_opt) {
                vector_memcpy(sol.x_opt, x);
                vector_memcpy(sol.v_opt, v);
            }
            sol.num_iters = iter + 1;
            break;
        }

        if (compute_newton_step(A, r_dual, r_primal, x, dx, dv, M, work, perm) != 0) {
            break;
        }

        backtracking_line_search(A, b, c, x, v, dx, dv, r_norm, x_new, v_new, r_dual, r_primal);

        vector_memcpy(x, x_new);
        vector_memcpy(v, v_new);
    }

    if (iter == MAX_ITER) {
        sol.status = FAILURE;
    }

    if (x) vector_free(x);
    if (v) vector_free(v);
    if (r_dual) vector_free(r_dual);
    if (r_primal) vector_free(r_primal);
    if (dx) vector_free(dx);
    if (dv) vector_free(dv);
    if (x_new) vector_free(x_new);
    if (v_new) vector_free(v_new);
    if (M) matrix_free(M);
    if (work) vector_free(work);
    if (perm) permutation_free(perm);

    return sol;
    return sol;
}

/**
 * @brief Solves a linear program with a feasible starting point using the interior point method.
 *
 * This function implements a path-following interior point method for solving a linear program of the form:
 *      minimize     cᵀx
 *      subject to   Ax = b
 *                   x > 0
 * starting from a strictly feasible point `x0`.
 *
 * It repeatedly solves a series of centering problems for increasing values of the barrier parameter `t`,
 * updating the solution estimate and tracking convergence via the duality gap.
 *
 * @param[in] A   Constraint matrix of size m × n.
 * @param[in] b   Right-hand side vector of size m.
 * @param[in] c   Cost vector of size n.
 * @param[in] x0  Strictly feasible initial point (x0 > 0), vector of size n.
 *
 * @return solution_t Structure containing:
 *         - status: OPTIMAL if convergence is reached, FAILURE otherwise.
 *         - x_opt: Optimal primal solution (if status == OPTIMAL).
 *         - v_opt: Optimal dual solution (if available).
 *         - opt_val: Optimal objective value cᵀx.
 *         - duality_gaps: Sequence of duality gap values for each iteration.
 *         - newton_steps: Number of Newton steps taken at each centering iteration.
 *
 * @note The caller is responsible for freeing the fields of the returned solution with `solution_free()`.
 * @warning If the centering step fails during any iteration, the solver terminates early with FAILURE status.
 */
static solution_t solve_feasible_start(const matrix_t *A, const vector_t *b, const vector_t *c, const vector_t *x0)
{
    LOG_INFO("solving feasible start");

    solution_t sol;
    vector_t *tc;
    vector_t *x;
    double gap;
    double t;

    solution_init(&sol);
    sol.duality_gaps = vector_alloc(MAX_ITER);
    sol.newton_steps = vector_alloc(MAX_ITER);
    tc = vector_alloc(c->size);
    x = vector_alloc(x0->size);
    vector_memcpy(x, x0);
    t = 1.0;

    for (sol.num_iters = 0; sol.num_iters < MAX_ITER; sol.num_iters++) {

        // Solve centering problem with current t
        solution_t center_sol;
        vector_memcpy(tc, c);
        vector_scale(tc, t);    // tc = t * c
        center_sol = solve_centering(A, b, tc, x);
        if (center_sol.status != OPTIMAL) {
            LOG_ERROR("Centering step failed at iteration %d", sol.num_iters);
            solution_free(&center_sol);
            break;
        }

        // Update solution
        vector_memcpy(x, center_sol.x_opt);

        // Compute and store duality gap and newton steps
        gap = A->size2 / t;
        vector_set(sol.duality_gaps, sol.num_iters, gap);
        vector_set(sol.newton_steps, sol.num_iters, center_sol.num_iters);

        // Check convergence
        if (gap < TOL) {
            sol.status = OPTIMAL;
            sol.x_opt = vector_alloc(x->size);
            sol.v_opt = vector_alloc(center_sol.v_opt->size);
            vector_memcpy(sol.x_opt, x);
            vector_memcpy(sol.v_opt, center_sol.v_opt);
            blas_ddot(c, x, &sol.opt_val);
            solution_free(&center_sol);
            break;
        }

        t = t * MU;
        solution_free(&center_sol);
    }

    LOG_INFO("Found %s solution", to_string(sol.status));

    vector_free(tc);
    vector_free(x);

    return sol;
}

/**
 * @brief Constructs the auxiliary LP problem for Phase I of the interior-point method.
 *
 * Given a linear system \( Ax = b \), this function constructs an auxiliary LP of the form:
 *
 * \[
 * \text{minimize} \quad t \\
 * \text{subject to} \quad [A, -A\mathbf{1}] \cdot [x; t] = b - A\mathbf{1}
 * \]
 * where \( \mathbf{1} \) is a vector of ones. This transformation allows solving an infeasible LP
 * by embedding it into a feasible auxiliary problem.
 *
 * @param[in]  A   Original constraint matrix of size m×n.
 * @param[in]  b   Original right-hand side vector of size m.
 * @param[out] A1  Pointer to the newly allocated auxiliary constraint matrix of size m×(n+1).
 * @param[out] b1  Pointer to the newly allocated auxiliary right-hand side vector of size m.
 * @param[out] c1  Pointer to the newly allocated auxiliary cost vector of size (n+1).
 *
 * @note The caller is responsible for freeing the output pointers (`*A1`, `*b1`, `*c1`)
 *       using the appropriate GSL functions.
 */
static void build_auxiliary_lp(const matrix_t *A, const vector_t *b, matrix_t **A1, vector_t **b1, vector_t **c1)
{
    LOG_INFO("Building Auxiliary LP");

    vector_t *A_ones;
    vector_t *ones;
    size_t rows;
    size_t cols;

    rows = A->size1;
    cols = A->size2;

    *A1 = matrix_alloc(rows, cols + 1);
    *b1 = vector_alloc(rows);
    *c1 = vector_alloc(cols + 1);
    A_ones = vector_alloc(rows);
    ones = vector_alloc(cols);

    // Build A1 = [A, -A*ones]
    vector_set_ones(ones);
    blas_dgemv(BLAS_NO_TRANSPOSE, -1.0, A, ones, 0.0, A_ones);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++)
            matrix_set(*A1, i, j, matrix_get(A, i, j));
        matrix_set(*A1, i, cols, vector_get(A_ones, i));
    }

    // Build b1 = b - A*ones
    vector_memcpy(*b1, b);
    blas_daxpy(-1.0, A_ones, *b1);

    // Build c1 = [0, ..., 0, 1]
    vector_set_zero(*c1);
    vector_set(*c1, cols, 1.0);

    if (A_ones) vector_free(A_ones);
    if (ones) vector_free(ones);
}

/**
 * @brief Builds the initial point for the auxiliary LP.
 *
 * Computes the initial feasible point \( z_0 \in \mathbb{R}^{n+1} \) for the auxiliary LP using:
 * \[
 * z_0 = [x + (t - 1) \cdot \mathbf{1}; \; t], \quad t = 2 - \min(x)
 * \]
 * This guarantees strict feasibility of the starting point in Phase I.
 *
 * @param[in]  x   The original (possibly infeasible) solution vector of size n.
 * @param[out] z0  Pointer to the newly allocated auxiliary starting point of size (n+1).
 *
 * @note The caller is responsible for freeing `*z0` using `vector_free()`.
 */
static void build_auxiliary_initial_point(const vector_t *x, vector_t **z0)
{
    LOG_INFO("Building Auxiliary initial point");

    vector_t *ones;
    vector_t *tmp;
    size_t v_len;
    double t;

    v_len = x->size;

    // Compute z0 = [x + (t-1)*ones, t], with t = 2 - min_x
    ones = vector_alloc(v_len);
    vector_set_ones(ones);
    tmp = vector_alloc(v_len);
    t = 2.0 - vector_min(x);
    vector_memcpy(tmp, x);
    blas_daxpy(t - 1.0, ones, tmp);
    *z0 = vector_concat(tmp, t);

    if (ones) vector_free(ones);
    if (tmp) vector_free(tmp);
}

/**
 * @brief Extracts a feasible starting point from the auxiliary LP optimal solution.
 *
 * Given an optimal solution \( z^* = [x^*; t^*] \in \mathbb{R}^{n+1} \) of the auxiliary LP,
 * this function returns the corrected feasible vector:
 * \[
 * x = x^* - (t^* - 1) \cdot \mathbf{1}
 * \]
 *
 * @param[in] aux_x_opt Optimal solution vector of the auxiliary LP (size n+1).
 * @param[in] v_len     Dimension of the original LP variables (i.e., n).
 *
 * @return A newly allocated GSL vector containing the feasible point for the original LP.
 *
 * @warning The caller is responsible for freeing the returned vector using `vector_free()`.
 */
static vector_t *get_auxiliary_opt(const vector_t *aux_x_opt, const size_t v_len)
{
    LOG_INFO("Extracting initial point from Auxiliary optimal solution");

    vector_t *x_opt;
    double t_opt;

    x_opt = vector_alloc(v_len);
    t_opt = vector_get(aux_x_opt, v_len);
    for (size_t i = 0; i < v_len; i++) {
        double val = vector_get(aux_x_opt, i) - (t_opt - 1.0);
        vector_set(x_opt, i, val);
    }

    return x_opt;
}

/**
 * @brief Solves the auxiliary LP problem to find a feasible starting point for the original LP.
 *
 * Constructs the auxiliary LP using `build_auxiliary_lp`, generates an initial feasible point with
 * `build_auxiliary_initial_point`, then solves the problem using `solve_feasible_start()`.
 * If the auxiliary LP is feasible and the optimal value is less than 1.0, a feasible point is extracted.
 *
 * @param[in] A Constraint matrix of the original LP (size m×n).
 * @param[in] b Right-hand side vector of the original LP (size m).
 * @param[in] x Infeasible solution estimate (size n), used to initialize the auxiliary LP.
 *
 * @return A `solution_t` structure containing:
 *         - status: `OPTIMAL` if a feasible point was found, `INFEASIBLE` otherwise.
 *         - x_opt:  Feasible starting point (if status is `OPTIMAL`).
 *         - Diagnostics (e.g., number of iterations, final value, etc.).
 *
 * @note The returned solution must be freed using `solution_free()` when no longer needed.
 */
static solution_t solve_auxiliary_lp(const matrix_t *A, const vector_t *b, const vector_t *x)
{
    LOG_INFO("Solving auxiliary problem");

    solution_t aux_sol;
    solution_t sol;
    matrix_t *A1;
    vector_t *b1;
    vector_t *c1;
    vector_t *z0;

    solution_init(&sol);
    build_auxiliary_lp(A, b, &A1, &b1, &c1);
    build_auxiliary_initial_point(x, &z0);

    aux_sol = solve_feasible_start(A1, b1, c1, z0);
    if (aux_sol.status == OPTIMAL && aux_sol.opt_val < 1.0) {
        LOG_INFO("Auxiliary Problem il feasible, extracting solution");
        sol.status = OPTIMAL;
        sol.x_opt = get_auxiliary_opt(aux_sol.x_opt, A->size2);
    } else {
        LOG_WARNING("Auxiliary Problem is INFEASIBLE - status=%d - opt_val=%f - n_iters=%d", sol.status, sol.opt_val, sol.num_iters);
        sol.status = INFEASIBLE;
    }

    if (A1) matrix_free(A1);
    if (b1) vector_free(b1);
    if (c1) vector_free(c1);
    if (z0) vector_free(z0);
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
 * @warning The caller is responsible for freeing the returned vector using `vector_free()`.
 * @note The pseudo-inverse is computed via SVD as \( A^+ = V \Sigma^+ U^T \) or \( A^+ = U \Sigma^+ V^T \)
 *       depending on the shape of A.
 */
static vector_t *pinv_solve(const matrix_t *A, const vector_t *b)
{
    LOG_INFO("Solving pseudo-inverse via SVD");

    vector_t *work;
    vector_t *tmp;
    vector_t *x;
    vector_t *S;
    matrix_t *U;
    matrix_t *V;
    size_t rows;
    size_t cols;

    rows = A->size1;
    cols = A->size2;

    if (rows >= cols) {
        U = matrix_alloc(rows, cols);
        V = matrix_alloc(cols, cols);
        work = vector_alloc(cols);
        tmp = vector_alloc(cols);
        S = vector_alloc(cols);
        x = vector_alloc(cols);

        matrix_memcpy(U, A);
        linalg_sv_decomp(U, V, S, work);

        blas_dgemv(BLAS_TRANSPOSE, 1.0, U, b, 0.0, tmp);    // tmp = Uᵗ b

        // Pseudoinverse: tmp = S⁺ * tmp
        for (size_t i = 0; i < cols; i++) {
            double s_val = vector_get(S, i);
            if (s_val > PINV_TOL) {
                vector_set(tmp, i, vector_get(tmp, i) / s_val);
            } else {
                vector_set(tmp, i, 0.0);
            }
        }

        blas_dgemv(BLAS_NO_TRANSPOSE, 1.0, V, tmp, 0.0, x);  // x = V * Σ⁺ * tmp
    } else {
        // Transpose A and compute SVD of Aᵗ
        U = matrix_alloc(cols, rows);
        V = matrix_alloc(rows, rows);
        work = vector_alloc(rows);
        tmp = vector_alloc(rows);
        S = vector_alloc(rows);
        x = vector_alloc(cols);

        matrix_transpose_memcpy(U, A);                  // U = Aᵗ
        linalg_sv_decomp(U, V, S, work);

        blas_dgemv(BLAS_TRANSPOSE, 1.0, V, b, 0.0, tmp);    // tmp = Vᵗ b

        // Pseudoinverse: tmp = S⁺ * tmp
        for (size_t i = 0; i < rows; i++) {
            double s_val = vector_get(S, i);
            if (s_val > PINV_TOL) {
                vector_set(tmp, i, vector_get(tmp, i) / s_val);
            } else {
                vector_set(tmp, i, 0.0);
            }
        }

        blas_dgemv(BLAS_NO_TRANSPOSE, 1.0, U, tmp, 0.0, x);  // x = V * Σ⁺ * tmp
    }

    if (work) vector_free(work);
    if (tmp) vector_free(tmp);
    if (U) matrix_free(U);
    if (V) matrix_free(V);
    if (S) vector_free(S);

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
static solution_t phase_one(const matrix_t *A, const vector_t *b)
{
    LOG_INFO("Solving phase I");

    solution_t sol;
    vector_t *x;

    solution_init(&sol);

    // Step 1: Compute initial point via pseudo-inverse
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
        sol.x_opt = vector_alloc(x->size);
        vector_memcpy(sol.x_opt, x);
    }

    if (x) vector_free(x);
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
solution_t solve(matrix_t *A, vector_t *b, vector_t *c)
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