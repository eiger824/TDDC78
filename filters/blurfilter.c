/*
File: blurfilter.c

Implementation of blurfilter function.

*/
#include <stdio.h>
#include "blurfilter.h"
#include "ppmio.h"


pixel* pix(pixel* image, const int xx, const int yy, const int xsize)
{
    register int off = xsize*yy + xx;

#ifdef DBG
    if (off >= MAX_PIXELS)
    {
        fprintf(stderr, "\n Terribly wrong: %d %d %d\n",xx,yy,xsize);
    }
#endif
    return (image + off);
}


#ifdef WITH_PTHREADS
void blurfilter_x(const int xstart, const int ystart,
        const int xend, const int yend, const int xsize,
        pixel* src, const int radius, const double *w)
{

    int x,y,x2,wi;
    double r,g,b,n,wc;

    for (y = ystart; y < yend; y++)
    {
        for (x = xstart; x < xend; x++)
        {
            r = w[0] * pix(src, x, y, xsize)->r;
            g = w[0] * pix(src, x, y, xsize)->g;
            b = w[0] * pix(src, x, y, xsize)->b;
            n = w[0];
            for ( wi = 1; wi <= radius; wi++)
            {
                wc = w[wi];
                x2 = x - wi;
                if (x2 >= 0)
                {
                    r += wc * pix(src, x2, y, xsize)->r;
                    g += wc * pix(src, x2, y, xsize)->g;
                    b += wc * pix(src, x2, y, xsize)->b;
                    n += wc;
                }
                x2 = x + wi;
                if (x2 < xend)
                {
                    r += wc * pix(src, x2, y, xsize)->r;
                    g += wc * pix(src, x2, y, xsize)->g;
                    b += wc * pix(src, x2, y, xsize)->b;
                    n += wc;
                }
            }
            pix(src,x,y, xsize)->r = r/n;
            pix(src,x,y, xsize)->g = g/n;
            pix(src,x,y, xsize)->b = b/n;
        }
    }
}

void blurfilter_y(const int xstart, const int ystart,
        const int xend, const int yend, const int xsize,
        pixel* src, const int radius, const double *w)
{

    int x,y,y2,wi;
    double r,g,b,n,wc;

    for (y=ystart; y < yend; y++)
    {
        for (x=xstart; x < xend; x++)
        {
            r = w[0] * pix(src, x, y, xsize)->r;
            g = w[0] * pix(src, x, y, xsize)->g;
            b = w[0] * pix(src, x, y, xsize)->b;
            n = w[0];
            for ( wi = 1; wi <= radius; wi++)
            {
                wc = w[wi];
                y2 = y - wi;
                if (y2 >= 0)
                {
                    r += wc * pix(src, x, y2, xsize)->r;
                    g += wc * pix(src, x, y2, xsize)->g;
                    b += wc * pix(src, x, y2, xsize)->b;
                    n += wc;
                }
                y2 = y + wi;
                if (y2 < yend)
                {
                    r += wc * pix(src, x, y2, xsize)->r;
                    g += wc * pix(src, x, y2, xsize)->g;
                    b += wc * pix(src, x, y2, xsize)->b;
                    n += wc;
                }
            }
            pix(src,x,y, xsize)->r = r/n;
            pix(src,x,y, xsize)->g = g/n;
            pix(src,x,y, xsize)->b = b/n;
        }
    }
}
#endif

#ifdef WITH_MPI

/* For the horizontal filter, every core would receive a portion of the image, so they can process it row by row */
void blurfilter_x(pixel * src, const uint nr_elems, const uint xsize, const uint radius, const double * w)
{
    int x,y,x2,y2,wi;
    double r,g,b,n, wc;
    uint i;

    for (i = 0; i< nr_elems; ++i)
    {
        y = i / xsize;
        x = i % xsize;

        r = w[0] * pix(src, x, y, xsize)->r;
        g = w[0] * pix(src, x, y, xsize)->g;
        b = w[0] * pix(src, x, y, xsize)->b;
        n = w[0];

        for ( wi = 1; wi <= radius; wi++)
        {
            wc = w[wi];
            x2 = x - wi;
            if (x2 >= 0)
            {
                r += wc * pix(src, x2, y, xsize)->r;
                g += wc * pix(src, x2, y, xsize)->g;
                b += wc * pix(src, x2, y, xsize)->b;
                n += wc;
            }
            x2 = x + wi;
            if (x2 < xsize)
            {
                r += wc * pix(src, x2, y, xsize)->r;
                g += wc * pix(src, x2, y, xsize)->g;
                b += wc * pix(src, x2, y, xsize)->b;
                n += wc;
            }
        }
        pix(src,x,y, xsize)->r = r/n;
        pix(src,x,y, xsize)->g = g/n;
        pix(src,x,y, xsize)->b = b/n;
    }
}
/* For the vertical filter, ROOT needs to send A COPY of the whole image to every core, and they will select
 * the relevant pixels to process according to the function below */
void blurfilter_y(const uint xstart, const uint ystart,
        const uint xend, const uint yend,
        constconst uint xsize, pixel * src, const uint radius, const double * w)
{
    int x,y,x2,y2,wi;
    double r,g,b,n, wc;
    uint i;

    for (x = xstart; x < xend; ++x)
    {
        for (y = ystart; y < yend; ++y)
        {

            r = w[0] * pix(src, x, y, xsize)->r;
            g = w[0] * pix(src, x, y, xsize)->g;
            b = w[0] * pix(src, x, y, xsize)->b;
            n = w[0];

            for ( wi = 1; wi <= radius; wi++)
            {
                wc = w[wi];
                x2 = x - wi;
                if (x2 >= 0)
                {
                    r += wc * pix(src, x2, y, xsize)->r;
                    g += wc * pix(src, x2, y, xsize)->g;
                    b += wc * pix(src, x2, y, xsize)->b;
                    n += wc;
                }
                x2 = x + wi;
                if (x2 < xsize)
                {
                    r += wc * pix(src, x2, y, xsize)->r;
                    g += wc * pix(src, x2, y, xsize)->g;
                    b += wc * pix(src, x2, y, xsize)->b;
                    n += wc;
                }
            }
            pix(src,x,y, xsize)->r = r/n;
            pix(src,x,y, xsize)->g = g/n;
            pix(src,x,y, xsize)->b = b/n;
        }
    }
}
#endif
