/*
  File: blurfilter.h

  Declaration of pixel structure and blurfilter function.
    
 */

#ifndef _BLURFILTER_H_
#define _BLURFILTER_H_

#include "defs.h"

void blurfilter_x(const int xstart, const int ystart,
        const int xend, const int yend, const int xsize,
        pixel* src, const int radius, const double *w);

void blurfilter_y(const int xstart, const int ystart,
        const int xend, const int yend, const int xsize,
        pixel* src, const int radius, const double *w);

#endif
