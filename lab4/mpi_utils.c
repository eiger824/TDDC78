#include <stdio.h>
#include <math.h>

#include "mpi_utils.h"
#include "definitions.h"

int compute_grid_dimensions(const uint nproc, int * dims)
{
    /* Superweird case */
    if (nproc == 0)
        return 1;

    uint i;
    if (check_if_prime(nproc)) 
    {
        *dims        = 1;        /* 1 single row */
        *(dims + 1)  = nproc;    /* As many columns as processors */
    }
    else
    {
        // Compute the highest factors
        compute_closest_factors(nproc, dims);
    }
    return 0;
}

int check_if_prime(const uint nproc)
{
    uint i;
    for (i = 2; i <= nproc / 2; ++i)
        if ( !(nproc % i) )
            return 0;
    return 1;
}

int compute_closest_factors(const uint nproc, int * dims)
{
    uint start = sqrt(nproc);
    while (nproc % start != 0)
        start--;
    *dims       = start;
    *(dims + 1) = nproc / start;
    return 0;
}

// Return coordinates of the grid region
cord_t * get_my_grid_boundaries(const uint hsize, const uint vsize, const int * grid, const int * dims)
{
    cord_t * out = (cord_t * ) malloc (sizeof(cord_t));

    out->x0 = grid[0] * hsize / dims[1];
    out->y0 = grid[1] * vsize / dims[0];

    out->x1 = (grid[0] + 1) * hsize / dims[1];;
    out->y1 = (grid[1] + 1) * vsize / dims[0];

    return out;
}

bool is_particle_in_grid_boundary(pcord_t * p1, cord_t * limits) 
{
    return p1->x == limits->x0 || p1->x == limits->x1 || p1->y == limits->y0 || p1->y == limits->y1; 
}

void print_limits(uint id, cord_t * limits)
{
    printf("[ID=%u] (x0,y0) = (%.2f,%.2f),\t\t(x1,y1) = (%.2f,%.2f)\n", id,
            limits->x0, limits->y0, limits->x1, limits->y1);
}
