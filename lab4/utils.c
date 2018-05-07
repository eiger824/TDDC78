#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/ioctl.h>

#include "utils.h"
#include "log.h"

static const int x_offset = 3;
static const int y_offset = 5;

void print_box(const uint horiz_size, const uint vert_size, pcord_t * particles, const uint nr_particles, const int * dims)
{
    struct winsize w;
    int err = ioctl(0, TIOCGWINSZ, &w);
    int rows;
    int cols;
    if (err == -1)
    {
        // In case of error, just set a 80x24 size (standard)
        rows = 24;
        cols = 80;
    }
    else
    {
        rows = (w.ws_row / 2) - y_offset;
        cols = w.ws_col - x_offset;
    }
    // Set the relative dimensions as well
    if (rows < cols)
    {
        cols = rows * horiz_size / vert_size; 
    }
    else if (cols < rows)
    {
        rows = cols * vert_size / horiz_size;
    }

    // Make cols and rows be divisible by the dimensions
    int ystep = cols / dims[1];
    if (ystep != 0)
        cols += (dims[1] - cols % dims[1]);
    int xstep = rows / dims[0];
    if (xstep != 0)
        rows += (dims[0] - rows % dims[0]);

    printf("Scaled box representation (cols=%d, rows=%d) (grid %dx%d):\n ", cols, rows, dims[0], dims[1]);

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
    // Print it
    for (i = 0; i < 2*cols; ++i)
        printf("~");
    printf(">X\n|");

    int toggle = 0;
    for (i = 0; i < rows; ++i)
    {
        for (j = 0; j < cols; ++j)
        {
            if (j % (cols/dims[1]) != 0)
                printf("%c%s", (!toggle ? box[i][j] : '~'), (j < cols - 1 ? " ": ""));
            else
                if (j > 0)
                    printf("%c%s", (!toggle ? box[i][j] : '~'), (j < cols - 1 ? "|": ""));
        }
        printf(" %s|\n", (!toggle ? ".": "~"));

        if (i > 0 && i % (rows/dims[0]) == 0)
            toggle = 1;
        else
            toggle = 0;
        // Beginning of a new row
        printf("|");
    }
    for (i = 0; i < 2*cols; ++i)
        printf("~");
    printf("\nY\n");

    // Free stuff
    for (i = 0; i < rows; ++i)
        free(*(box + i));
    free(box);
}
