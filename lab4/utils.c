#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/ioctl.h>

#include "utils.h"
#include "log.h"

static const int x_offset = 3;
static const int y_offset = 5;

void print_box(const uint horiz_size, const uint vert_size, pcord_t * particles, const uint nr_particles)
{
    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);
    int rows = (w.ws_row / 2) - y_offset;
    int cols = w.ws_col - x_offset;
    // Set the relative dimensions as well
    if (rows < cols)
    {
        cols = rows * horiz_size / vert_size; 
    }
    else if (cols < rows)
    {
        rows = cols * vert_size / horiz_size;
    }
    log_debug("Printing a %dx%d box", rows, cols);

    // Init the box 
    uint i,j;
    char ** box = (char **) malloc (sizeof (char*) * rows);
    for (i = 0; i < rows; ++i)
    {
        *(box + i) = (char *) malloc(sizeof (char) * cols);
        memset(*(box+i), '.', cols);
    }
    // Fill the scaled representation
    for (i = 0; i < nr_particles; ++i)
    {
        int rela_x = particles[i].x * cols / horiz_size;
        int rela_y = particles[i].y * rows / vert_size;
        box[rela_y][rela_x] = 'o';
    }
    printf("Scaled box representation:\n ");
    // Print it
    for (i = 0; i < cols; ++i)
        printf("~");
    printf(">X\n|");

    for (i = 0; i < rows; ++i)
    {
        for (j = 0; j < cols; ++j)
        {
            printf("%c", box[i][j]);
        }
        printf("|\n|");
    }
    for (i = 0; i < cols; ++i)
        printf("~");
    printf("\nY\n");

    // Free stuff
    for (i = 0; i < rows; ++i)
        free(*(box + i));
    free(box);
}
