/*
  File: thresfilter.h

  Declaration of pixel structure and thresfilter function.
    
 */
#ifndef _THRESFILTER_H_
#define _THRESFILTER_H_

/* NOTE: This structure must not be padded! */
typedef struct _pixel
{
    unsigned char r,g,b;
} pixel;

typedef unsigned int uint;

uint get_px_sum(pixel * src, const uint nr_elems);
void thresfilter(pixel* src, const uint nr_elems, const uint average);


#endif
