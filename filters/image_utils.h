/*
 * =====================================================================================
 *
 *       Filename:  image_utils.h
 *
 *    Description:  Utilities to work with images
 *
 *        Version:  1.0
 *        Created:  2018-04-23 11:57:29
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Santiago Pagola (), santipagola@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef IMG_UTILS_H_
#define IMG_UTILS_H_

#include "defs.h"

struct matrix
{
    uint xsize;
    uint ysize;
    pixel ** m;
};

struct matrix * from_array(pixel * src, const uint xsize, const uint ysize, struct matrix * mt);
pixel * to_array(struct matrix * mt, pixel * out);
void print_matrix(struct matrix * mt);
void print_array(pixel * src, const uint nr_elems);
void init_matrix(struct matrix * mt);
void fill_matrix(struct matrix * mt);
struct matrix * adj_matrix(struct matrix * orig, struct matrix * adj);
void free_matrix(struct matrix * mt);

#endif /* IMG_UTILS_H_ */
