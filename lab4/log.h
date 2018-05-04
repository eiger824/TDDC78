#ifndef LOG_H_
#define LOG_H_

#include <stdio.h>
#include <stdarg.h>
#include <stdbool.h>

#define ERROR               0
#define INFO                1

#define log_info(...)       print(INFO, __FILE__, __LINE__, __VA_ARGS__)
#define log_error(...)      print(ERROR, __FILE__, __LINE__, __VA_ARGS__)
#define log_set(...)        log_enable(__VA_ARGS__)

#define ERROR_COLOR         "\x1B[31m"
#define NORMAL_COLOR        "\x1B[0m"

void log_enable(bool enable);  /* Enables logging or not */

void print(unsigned level,       /* level: ERROR (stderr) / INFO (stdout) */
           const char* filename, /* filename: The current file where the logging line is */
           unsigned line,        /* line: The current line where the logging line is */
           const char* msg, ...); /* msg: Formated message with following (variadic) args */

#endif /*LOG_H_*/

