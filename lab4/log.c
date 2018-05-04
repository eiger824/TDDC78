#include <string.h>
#include <time.h>
#include <sys/time.h>

#include "log.h"

static struct   timeval time_before, time_after, time_result;
static bool     g_enabled = false;

void log_enable(bool enable)
{
    g_enabled = enable;
}

void print(unsigned level,       /* level: ERROR (stderr) / INFO (stdout) */
        const char* filename, /* filename: The current file where the logging line is */
        unsigned line,        /* line: The current line where the logging line is */
        const char* msg, ...) /* msg: Formated message with following (variadic) args */
{
    if ( (g_enabled && level) || !level )  /* Log to stdout if verbose enabled ** or ** if logging to stderr ALWAYS */
    {
        //update the time struct
        gettimeofday(&time_after, NULL);

        time_result.tv_sec = time_after.tv_sec - time_before.tv_sec;
        time_result.tv_usec = time_after.tv_usec - time_before.tv_usec;

        va_list args;

        long int usecs = (long int)time_result.tv_usec;
        time_t timer;
        struct tm* tm_info;

        time(&timer);
        tm_info = localtime(&timer);

        char time_str[50];
        time_str[strftime(time_str, 26, "%Y-%m-%d %H:%M:%S", tm_info)] = '\0';

        char usecs_str[20];
        usecs_str[sprintf(usecs_str, "%06ld", usecs)] = '\0';

        if (!level)
            fprintf(stderr, ERROR_COLOR"[%s.%s] (%s:%d) ", time_str, usecs_str, filename, line);
        else
            printf("[%s%s] (%s:%d) ", time_str, usecs_str, filename, line);

        va_start(args, msg);
        if (!level)
        {
            vfprintf(stderr, msg, args);
            fprintf(stderr, NORMAL_COLOR"\n");
        }
        else
        {
            vfprintf(stdout, msg, args);
            printf(NORMAL_COLOR"\n");
        }
        va_end(args);

        //and update time struct
        memcpy(&time_before, &time_after, sizeof(time_before));
    }
}
