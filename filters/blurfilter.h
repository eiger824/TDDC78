/*
  File: blurfilter.h

  Declaration of pixel structure and blurfilter function.
    
 */

#ifndef _BLURFILTER_H_
#define _BLURFILTER_H_

typedef unsigned int uint;

/* NOTE: This structure must not be padded! */
typedef struct _pixel {
    unsigned char r,g,b;
} pixel;

#ifdef WITH_PTHREADS
void blurfilter_x(const int xstart, const int ystart,
        const int xend, const int yend, const int xsize,
        pixel* src, const int radius, const double *w);

void blurfilter_y(const int xstart, const int ystart,
        const int xend, const int yend, const int xsize,
        pixel* src, const int radius, const double *w);
#endif

#ifdef WITH_MPI
void blurfilter_x(pixel * src, const uint nr_elems, const uint xsize, const uint radius, const double * w);
void blurfilter_y(const uint xstart, const uint ystart,
        const uint xend, const uint yend,
        const uint xsize, pixel * src, const uint radius, const double * w);
#endif

#endif
