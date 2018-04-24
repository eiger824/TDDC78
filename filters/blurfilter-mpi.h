/*
  File: blurfilter.h

  Declaration of pixel structure and blurfilter function.
    
 */

#ifndef _BLURFILTER_H_
#define _BLURFILTER_H_

#include "defs.h"

void blurfilter(pixel * src, const uint nr_elems, const uint xsize, const uint radius, const double * w);

#endif
