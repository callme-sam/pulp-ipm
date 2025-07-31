#ifndef UTILS_H_
#define UTILS_H_

#include <stdbool.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef enum {
    LOG_LEVEL_QUIET   = 0,
    LOG_LEVEL_ERROR   = 1,
    LOG_LEVEL_WARNING = 2,
    LOG_LEVEL_INFO    = 3,
    LOG_LEVEL_DEBUG   = 4
} LogLevel;

#ifndef CURRENT_LOG_LEVEL
#define CURRENT_LOG_LEVEL LOG_LEVEL_INFO
#endif

#define LOG(level, fmt, ...) do { \
    if ((level) <= CURRENT_LOG_LEVEL) { \
        const char* level_str; \
        FILE* stream = stdout; \
        switch(level) { \
            case LOG_LEVEL_ERROR:   level_str = "ERROR";   stream = stderr; break; \
            case LOG_LEVEL_WARNING: level_str = "WARNING"; stream = stderr; break; \
            case LOG_LEVEL_INFO:    level_str = "INFO";    break; \
            case LOG_LEVEL_DEBUG:   level_str = "DEBUG";   break; \
            default:               level_str = "UNKNOWN"; break; \
        } \
        fprintf(stream, "[%s] %s:%d - " fmt "\n", \
                level_str, __func__, __LINE__, ##__VA_ARGS__); \
    } \
} while(0)

#define LOG_ERROR(fmt, ...)   LOG(LOG_LEVEL_ERROR,   fmt, ##__VA_ARGS__)
#define LOG_WARNING(fmt, ...) LOG(LOG_LEVEL_WARNING, fmt, ##__VA_ARGS__)
#define LOG_INFO(fmt, ...)    LOG(LOG_LEVEL_INFO,    fmt, ##__VA_ARGS__)
#define LOG_DEBUG(fmt, ...)   LOG(LOG_LEVEL_DEBUG,   fmt, ##__VA_ARGS__)

void fill_matrix(gsl_matrix *M);
void fill_vector(gsl_vector *v);

void print_matrix(const gsl_matrix *A);
void print_vector(const gsl_vector *v);

void log_matrix(LogLevel level, const gsl_matrix *m, const char *name);
void log_vector(LogLevel level, const gsl_vector *v, const char *name);

bool has_nonpositive_elements(const gsl_vector *v);

gsl_vector *vector_ones(size_t n);
gsl_vector *vector_concat(const gsl_vector *v, double b);

void generate_lp(gsl_matrix **A, gsl_vector **b, gsl_vector **c);

#endif  /* UTILS_H_ */