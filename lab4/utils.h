#ifndef UTILS_H_
#define UTILS_H_

#include "definitions.h"

void print_box(const uint horiz_size, const uint vert_size, pcord_t * particles, const uint nr_particles, const int * dims);

void set_bogus_values(pcord_t ** matrix, int rows, int cols);

#endif /* UTILS_H_ */
