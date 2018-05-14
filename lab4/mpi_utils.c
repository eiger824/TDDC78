#include <stdio.h>
#include <math.h>

#include "mpi_utils.h"
#include "definitions.h"
#include "log.h"

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

    out->x0 = grid[1] * hsize / dims[1];
    out->y0 = grid[0] * vsize / dims[0];

    out->x1 = (grid[1] + 1) * hsize / dims[1];
    out->y1 = (grid[0] + 1) * vsize / dims[0];

    return out;
}

bool is_particle_inside_grid_boundary(pcord_t * p1, cord_t * my_limits) 
{
    return p1->x >= my_limits->x0 && p1->x <= my_limits->x1 && p1->y >= my_limits->y0 && p1->y <= my_limits->y1; 
}

bool is_particle_outside_grid_boundary(pcord_t * p1, cord_t * my_limits, nloc_t * dir)
{
    if (!is_particle_inside_grid_boundary(p1, my_limits))
    {
        // Update the neighbor where it actually is
        get_particle_dir(p1, my_limits, dir);
        return true;
    }
    return false;
}

void print_limits(uint id, cord_t * my_limits)
{
    log_warning("[ID=%02u] (x0,y0) = (%.2f,%.2f),\t(x1,y1) = (%.2f,%.2f)", id,
            my_limits->x0, my_limits->y0, my_limits->x1, my_limits->y1);
}

// Given a coordinate p1, determine its cartesian coordinates in the MPI 2D topology
void get_particle_dir(pcord_t * p, cord_t * my_limits, nloc_t * dir)
{
    int x = p->x;
    int y = p->y;

    if (x < my_limits->x0) *dir = LEFT;
    else if (x > my_limits->x1) *dir = RIGHT;
    else if (y < my_limits->y0) *dir = TOP;
    else if (y > my_limits->y1) *dir = BOTTOM;
}
