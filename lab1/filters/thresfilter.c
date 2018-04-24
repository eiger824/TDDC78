#include "thresfilter.h"

uint get_px_sum(pixel * src, const uint nr_elems)
{
    uint sum = 0;
    uint i;

    for (i = 0; i < nr_elems; i++)
    {
        sum += (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
    }

    return sum;
}

void thresfilter(pixel* src, const uint nr_elems, const uint average)
{
    uint i, psum;

    for (i = 0; i < nr_elems; i++)
    {
        psum = (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
        if( psum < average)
        {
            src[i].r = src[i].g = src[i].b = 0;
        }
        else
        {
            src[i].r = src[i].g = src[i].b = 255;
        }
    }
}
