#ifndef LOG_H_
#define LOG_H_

#include <stdio.h>
#include <stdarg.h>
#include <stdbool.h>

typedef enum log_level
{
    SILENT = -1,
    ERROR,
    WARNING,
    INFO,
    DEBUG
} level_t;

#define log_debug(...)      print(DEBUG, __FILE__, __LINE__, __VA_ARGS__)
#define log_warning(...)    print(WARNING, __FILE__, __LINE__, __VA_ARGS__)
#define log_info(...)       print(INFO, __FILE__, __LINE__, __VA_ARGS__)
#define log_error(...)      print(ERROR, __FILE__, __LINE__, __VA_ARGS__)
#define log_set(...)        log_enable(__VA_ARGS__)

#define ERROR_COLOR         "\x1B[31m"
#define NORMAL_COLOR        "\x1B[0m"

void log_enable(level_t level);  /* Enables logging or not */

void print(level_t level,       /* level: ERROR (stderr) / INFO (stdout) */
           const char* filename, /* filename: The current file where the logging line is */
           unsigned line,        /* line: The current line where the logging line is */
           const char* msg, ...); /* msg: Formated message with following (variadic) args */

#endif /*LOG_H_*/

