#include "log.h"

/**
 * @brief Logs the contents of a GSL matrix with a custom label.
 *
 * This function prints the contents of a GSL matrix to the log output,
 * formatted as a human-readable 2D array with a descriptive label.
 * Each row is printed on a separate line, with elements separated by commas:
 * ```
 * Matrix <name> (size = rows x cols):
 *     [a_00, a_01, ..., a_0n]
 *     [a_10, a_11, ..., a_1n]
 *     ...
 * ```
 * Logging only occurs if the specified `level` is less than or equal to `CURRENT_LOG_LEVEL`.
 *
 * @param[in] level Logging severity level (e.g., DEBUG, INFO, ERROR).
 * @param[in] m Pointer to the GSL matrix to be logged.
 * @param[in] name Descriptive name of the matrix, shown in the log output.
 *
 * @note The function dynamically allocates memory to format the matrix string.
 *       It automatically releases the memory before returning.
 *
 * @warning If memory allocation fails, the function logs an error and returns without printing.
 */
void log_matrix(LogLevel level, const matrix_t *m, const char *name) {
    if (level > CURRENT_LOG_LEVEL) return;

    const size_t per_elem_len = 24;
    const size_t header_len = 64;
    const size_t indent_len = 4;
    size_t row_len = indent_len + per_elem_len * m->size2 + 2;
    size_t buf_size = header_len + row_len * m->size1 + 1;

    char *buf = malloc(buf_size);
    if (!buf) {
        LOG_ERROR("Failed to allocate memory for matrix logging");
        return;
    }

    size_t offset = snprintf(buf, buf_size, "Matrix %s (size = %zux%zu):\n", name, m->size1, m->size2);

    for (size_t i = 0; i < m->size1; ++i) {
        offset += snprintf(buf + offset, buf_size - offset, "\t[");
        for (size_t j = 0; j < m->size2; ++j) {
            offset += snprintf(buf + offset, buf_size - offset, "%g%s",
                               matrix_get(m, i, j), (j < m->size2 - 1) ? ", " : "");
        }
        offset += snprintf(buf + offset, buf_size - offset, "]\n");
    }
    buf[offset] = '\0';

    LOG(level, "%s", buf);
    free(buf);
}

/**
 * @brief Logs the contents of a GSL vector with a custom label.
 *
 * This function formats and prints the elements of a GSL vector using a specified log level.
 * The vector is printed in the format:
 * ```
 * Vector <name> (size = n):
 *     [e0, e1, ..., en]
 * ```
 * Logging only occurs if the specified level is less than or equal to the current log level.
 *
 * @param[in] level Logging severity level
 * @param[in] v Pointer to the GSL vector to be printed.
 * @param[in] name Descriptive name to be shown in the log output.
 *
 * @note The function allocates temporary memory to construct the formatted log message.
 *       It automatically frees the memory before returning.
 *
 * @warning If memory allocation fails, the function logs an error and exits early.
 */
void log_vector(LogLevel level, const vector_t *v, const char *name) {
   if (level > CURRENT_LOG_LEVEL) return;

   const size_t per_elem_len = 24;
    const size_t header_len = 64;
    size_t buf_size = header_len + per_elem_len * v->size;
    char *buf = malloc(buf_size);
    if (!buf) {
        LOG_ERROR("Failed to allocate memory for vector logging");
        return;
    }

    size_t offset = snprintf(buf, buf_size, "Vector %s (size = %zu):\n\t[", name, v->size);

    for (size_t i = 0; i < v->size; ++i) {
        offset += snprintf(buf + offset, buf_size - offset, "%g%s",
                           vector_get(v, i), (i < v->size - 1) ? ", " : "");
    }

    snprintf(buf + offset, buf_size - offset, "]");

    LOG(level, "%s", buf);

    free(buf);
}