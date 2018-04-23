/*
 * =====================================================================================
 *
 *       Filename:  defs.h
 *
 *    Description:  Common definitions shared across translation units.
 *
 *        Version:  1.0
 *        Created:  2018-04-23 13:10:35
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Santiago Pagola (), santipagola@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef DEFS_H_
#define DEFS_H_

typedef unsigned int uint;

/* NOTE: This structure must not be padded! */
typedef struct _pixel
{
    unsigned r,g,b;
} pixel;

#endif
