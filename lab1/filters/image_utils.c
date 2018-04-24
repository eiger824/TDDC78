/*
 * =====================================================================================
 *
 *       Filename:  image_utils.c
 *
 *    Description:  Utilities to work with images - implementation file
 *
 *        Version:  1.0
 *        Created:  2018-04-23 12:01:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Santiago Pagola (), santipagola@gmail.com
 *   Organization:
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "image_utils.h"

void print_array(pixel * src, const uint nr_elems)
{
    uint i;
    pixel p;
    printf("{ ");
    for (i = 0; i < nr_elems; ++i)
    {
        p = *(src + i);
        printf("[%d,%d,%d] %s", p.r, p.g, p.b,
                (i < nr_elems - 1 ) ? ",":"");
    }
    printf("}\n");
}

void print_matrix(struct matrix * mt)
{
    uint i,x,y;
    pixel p;

    printf("\t  ");
    for (i = 0; i < mt->xsize; ++i)
        printf("%u\t\t", i);
    printf("\n\t  ");
    for (i = 0; i < mt->xsize; ++i)
        printf("--\t\t");

    printf("\n");
    for (y = 0; y < mt->ysize; ++y)
    {
        printf("%u\t|", y);
        for (x = 0; x < mt->xsize; ++x)
        {
            p = *(*(mt->m + y) + x);
            printf("%1s[%d,%d,%d]\t", " ", p.r, p.g, p.b);
        }
        printf("\n\n");
    }
}

void init_matrix(struct matrix * mt)
{
    uint i;

    mt->m = (pixel **) malloc(sizeof(pixel*) * mt->ysize);
    for (i = 0; i < mt->ysize; ++i)
    {
        *(mt->m + i) = (pixel *) malloc(sizeof(pixel) * mt->xsize);
    }
}

void fill_matrix(struct matrix * mt)
{
    uint cnt = 0;
    uint y,x;
    pixel p;

    for (y = 0; y < mt->ysize; ++y)
    {
        for (x = 0; x < mt->xsize; ++x)
        {
            p = *(*(mt->m + y) + x);
            p.r = ++cnt;
            p.g = cnt;
            p.b = cnt;
            *(*(mt->m + y) + x) = p;
        }
    }
}

/* WARNING: `src` must contain exactly mt->xsize * mt->ysize elements */
static void fill_matrix_from_array(struct matrix * mt, pixel * src, const uint nr_elems, const uint xsize)
{
    uint x,y;
    uint i;

    for (i = 0; i < nr_elems; ++i)
    {
        x = i % xsize;
        y = i / xsize;
        // Get the position (x,y) of the matrix and fill what's inside of src[i]
        *(*(mt->m + y) + x) = *(src + i);
    }
}


struct matrix * adj_matrix(struct matrix * orig, struct matrix * adj)
{
    uint x,y;

    /* Swap sizes */
    adj->xsize = orig->ysize;
    adj->ysize = orig->xsize;
    /* Init the adjoint matrix */
    init_matrix(adj);
    /* Swap the colums and rows */
    for (y = 0; y < adj->ysize; ++y)
    {
        for (x = 0; x < adj->xsize; ++x)
        {
            *(*(adj->m + y) + x) = *(*(orig->m + x) + y);
        }
    }
    /* For commodity purposes, return `adj` as well */
    return adj;
}

struct matrix * from_array(pixel * src, const uint xsize, const uint ysize, struct matrix * mt)
{
    mt->xsize = xsize;
    mt->ysize = ysize;
    init_matrix(mt);

    fill_matrix_from_array(mt, src, xsize * ysize, xsize);

    return mt;
}

pixel * to_array(struct matrix * mt, pixel * out)
{
    uint i;

    for (i = 0; i < mt->ysize; ++i)
    {
        memcpy(out + i * mt->xsize ,
                mt->m[i],
                mt->xsize * sizeof(pixel));
    }
    return out;
}

void free_matrix(struct matrix * mt)
{
    uint i;

    for (i = 0; i < mt->ysize; ++i)
        free(*(mt->m + i));

    free(mt->m);
}

