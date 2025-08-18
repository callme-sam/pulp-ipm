#include <stdbool.h>
#include <stdint.h>

#include "log.h"

static size_t utoa(uint64_t val, char *buf) {
    char temp[20];
    int i = 0;
    do {
        temp[i++] = '0' + val % 10;
        val /= 10;
    } while (val);
    for (int j = 0; j < i; ++j) {
        buf[j] = temp[i - j - 1];
    }
    return i;
}

static size_t ftoa(double val, char *buf) {
    char *start = buf;
    if (val < 0) {
        *buf++ = '-';
        val = -val;
    }
    uint64_t int_part = (uint64_t)val;
    double frac = val - int_part;

    buf += utoa(int_part, buf);
    *buf++ = '.';

    frac *= 1000000;
    uint64_t frac_part = (uint64_t)(frac + 0.5);

    uint64_t pow10 = 100000;
    while (frac_part < pow10 && pow10 > 1) {
        *buf++ = '0';
        pow10 /= 10;
    }

    buf += utoa(frac_part, buf);
    return buf - start;
}

/**
 * @brief Logs the contents of a matrix with a custom label.
 *
 * This function prints the contents of a matrix to the log output,
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
 * @param[in] m Pointer to the matrix to be logged.
 * @param[in] name Descriptive name of the matrix, shown in the log output.
 *
 * @note The function dynamically allocates memory to format the matrix string.
 *       It automatically releases the memory before returning.
 *
 * @warning If memory allocation fails, the function logs an error and returns without printing.
 */
void log_matrix(LogLevel level, const matrix_t *m, const char *name)
{
    if (level > CURRENT_LOG_LEVEL) return;

    const size_t per_elem_len = 32;
    const size_t header_len = 64;
    const size_t indent_len = 4;
    const size_t row_len = indent_len + per_elem_len * m->size2 + 2;
    const size_t buf_size = header_len + row_len * m->size1 + 1;

    char *buf = malloc(buf_size);
    if (!buf) {
        LOG_ERROR("Failed to allocate memory for matrix logging");
        return;
    }

    char *p = buf;

    const char *prefix = "Matrix ";
    const char *midfix = " (size = ";
    const char *xfix = "x";
    const char *suffix = "):\n";

    for (const char *s = prefix; *s; ++s) *p++ = *s;
    for (const char *s = name; *s; ++s) *p++ = *s;
    for (const char *s = midfix; *s; ++s) *p++ = *s;

    p += utoa(m->size1, p);
    for (const char *s = xfix; *s; ++s) *p++ = *s;
    p += utoa(m->size2, p);
    for (const char *s = suffix; *s; ++s) *p++ = *s;

    for (size_t i = 0; i < m->size1; ++i) {
        *p++ = '\t';
        *p++ = '[';
        for (size_t j = 0; j < m->size2; ++j) {
            p += ftoa(matrix_get(m, i, j), p);
            if (j < m->size2 - 1) {
                *p++ = ',';
                *p++ = ' ';
            }
        }
        *p++ = ']';
        *p++ = '\n';
    }

    *p = '\0';
    LOG(level, "%s", buf);
    free(buf);
}

/**
 * @brief Logs the contents of a vector with a custom label.
 *
 * This function formats and prints the elements of a vector using a specified log level.
 * The vector is printed in the format:
 * ```
 * Vector <name> (size = n):
 *     [e0, e1, ..., en]
 * ```
 * Logging only occurs if the specified level is less than or equal to the current log level.
 *
 * @param[in] level Logging severity level
 * @param[in] v Pointer to the vector to be printed.
 * @param[in] name Descriptive name to be shown in the log output.
 *
 * @note The function allocates temporary memory to construct the formatted log message.
 *       It automatically frees the memory before returning.
 *
 * @warning If memory allocation fails, the function logs an error and exits early.
 */
void log_vector(LogLevel level, const vector_t *v, const char *name) {
    if (level > CURRENT_LOG_LEVEL) return;

    const size_t per_elem_len = 32;
    const size_t header_len = 64;
    const size_t buf_size = header_len + per_elem_len * v->size + 3; // +3 for "[]\0"

    char *buf = malloc(buf_size);
    if (!buf) {
        LOG_ERROR("Failed to allocate memory for vector logging");
        return;
    }

    char *p = buf;

    // Header: "Vector <name> (size = N):\n\t["
    const char *prefix = "Vector ";
    const char *midfix = " (size = ";
    const char *suffix = "):\n\t[";

    for (const char *s = prefix; *s; ++s) *p++ = *s;
    for (const char *s = name; *s; ++s) *p++ = *s;
    for (const char *s = midfix; *s; ++s) *p++ = *s;
    p += utoa(v->size, p);
    for (const char *s = suffix; *s; ++s) *p++ = *s;

    for (size_t i = 0; i < v->size; ++i) {
        p += ftoa(vector_get(v, i), p);
        if (i < v->size - 1) {
            *p++ = ',';
            *p++ = ' ';
        }
    }

    *p++ = ']';
    *p = '\0';

    LOG(level, "%s", buf);
    free(buf);
}