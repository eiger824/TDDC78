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

