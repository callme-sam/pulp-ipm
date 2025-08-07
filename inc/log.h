#ifndef LOG_H_
#define LOG_H_

#include <stdio.h>

#include "matrix.h"
#include "vector.h"

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

void log_matrix(LogLevel level, const matrix_t *m, const char *name);
void log_vector(LogLevel level, const vector_t *v, const char *name);

#endif  /* LOG_H_ */