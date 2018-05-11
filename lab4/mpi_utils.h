#ifndef MPI_UTILS_H_
#define MPI_UTILS_H_

#include <stdarg.h>
#include <stdbool.h>

#include "definitions.h"

/*
 * Declare convenience macros defining the custom mpi datatype
 * needed for this lab.
 */


/*
 * DECLARE_MPI_COORDINATE_DATATYPE
 * Declares an MPI Datatype of the type "cord_t"
 *
 * X   :   Original type that we want to convert.
 * Y   :   MPI_Datatype, i.e., conversion of X.
 *
 */

#define     DECLARE_MPI_COORDINATE_DATATYPE(X,Y) \
    X item;\
MPI_Datatype Y;\
\
int block_lengths [] = {1, 1, 1, 1};\
MPI_Datatype block_types [] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};\
MPI_Aint start, displ[4];\
\
MPI_Address( &item, &start );\
MPI_Address( &item.x0, &displ[0] );\
MPI_Address( &item.x1, &displ[1] );\
MPI_Address( &item.y0, &displ[2] );\
MPI_Address( &item.y1, &displ[3] );\
\
displ[0] -= start;\
displ[1] -= start;\
displ[2] -= start;\
displ[3] -= start;\
\
MPI_Type_struct( 4, block_lengths, displ, block_types, &Y );\
MPI_Type_commit( &Y );\


/*
 * DECLARE_MPI_PARTICLE_COORDINATE_DATATYPE
 * Declares an MPI Datatype of the type "pcord_t"
 *
 * X   :   Original type that we want to convert.
 * Y   :   MPI_Datatype, i.e., conversion of X.
 *
 */
#define     DECLARE_MPI_PARTICLE_COORDINATE_DATATYPE(X,Y) \
    X item_pc;\
MPI_Datatype Y;\
\
MPI_Address( &item_pc, &start );\
MPI_Address( &item_pc.x, &displ[0] );\
MPI_Address( &item_pc.y, &displ[1] );\
MPI_Address( &item_pc.vx, &displ[2] );\
MPI_Address( &item_pc.vy, &displ[3] );\
\
displ[0] -= start;\
displ[1] -= start;\
displ[2] -= start;\
displ[3] -= start;\
\
MPI_Type_struct( 4, block_lengths, displ, block_types, &Y );\
MPI_Type_commit( &Y );\


/*
 * DECLARE_MPI_PARTICLE_DATATYPE
 * Declares an MPI Datatype of the type "particle_t"
 *
 * X    :   Original type that we want to convert.
 * Y    :   MPI_Datatype, i.e., conversion of X.
 * Z    :   MPI_Datatype for the Coordinate structure, needed.
 *
 */
#define     DECLARE_MPI_PARTICLE_DATATYPE(X, Y, Z) \
    X item_p;\
MPI_Datatype Y;\
\
int block_lengths_p [] = {1, 1};\
MPI_Datatype block_types_p [] = {Z, MPI_INT};\
MPI_Aint start_p, displ_p[2];\
\
MPI_Address( &item_p, &start_p );\
MPI_Address( &item_p.pcord, &displ_p[0] );\
MPI_Address( &item_p.ptype, &displ_p[1] );\
\
displ_p[0] -= start_p;\
displ_p[1] -= start_p;\
\
MPI_Type_struct( 2, block_lengths_p, displ_p, block_types_p, &Y );\
MPI_Type_commit( &Y );\


/*
 * Function:	compute_grid_dimensions
 * Brief:	    Given a number of processors, it computes the size of the grid each processor will be assigned
 * @param nproc:Number of available processors
 * @param dims:	The dimensions array
 * Returns:	    0 on success
 */
int compute_grid_dimensions(const uint nproc, int * dims);

/*
 * Function:	check_if_prime
 * Brief:	    Given a positive integer, it determines whether the number is prime or not
 * @param nproc:The positive integer to check
 * Returns:	    1 if prime, 0 otherwise
 */
int check_if_prime(const uint nproc);

/*
 * Function:	compute_closest_factors
 * Brief:       Given a positive integer, it computes its closest factors
 * @param nproc:The integer to check
 * @param dims:	The array in which to store the factors
 * Returns:	    0 on success, 1 otherwise
 */
int compute_closest_factors(const uint nproc, int * dims);

/*
 * Function:	get_my_grid_boundaries
 * Brief:	    Given a some MPI cartesian coordinates, it returns its physical limits
 * @param hsize:Horizontal size of the box
 * @param vsize:Vertical size of the box
 * @param grid:	My cartesian coordinates in the MPI topology
 * @param dims:	The dimension of the MPI 2D topology matrix.
 * Returns:	    The limits of the grid region
*/
cord_t * get_my_grid_boundaries(const uint hsize, const uint vsize, const int * grid, const int * dims);


bool is_particle_inside_grid_boundary(pcord_t * p1, cord_t * my_limits);
bool is_particle_outside_grid_boundary(pcord_t * p1, cord_t * my_limits,
        const uint hsize, const uint vsize, const int * dims, int * nbr_coord );

void print_limits(uint id, cord_t * limits);

void get_grid_region_of_particle(pcord_t * p, const uint hsize, const uint ysize, const int * dims, int * my_grid );
#endif  /* MPI_UTILS_H_ */

