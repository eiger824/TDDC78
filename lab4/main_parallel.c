#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <getopt.h>

#include <mpi.h>

#include "coordinate.h"
#include "definitions.h"
#include "physics.h"
#include "log.h"
#include "mpi_utils.h"


//Feel free to change this program to facilitate parallelization.

void help(const char * program)
{
    printf("Usage: %s [h|n|x|y] <sim-time>\n", program);
    printf("-h\tPrint this help and exit\n");
    printf("-n nr\tSet number of particles each processor will handle\n\t(default: %d, max: %d)\n",
            INIT_NO_PARTICLES, MAX_NO_PARTICLES);
    printf("-x size\tSet box horizontal size (default: %.2f)\n",
            BOX_HORIZ_SIZE);
    printf("-y size\tSet box vertical size (default: %.2f)\n",
            BOX_VERT_SIZE);
    printf("-v\tSet verbose logging\n");
}

float rand1()
{
    return (float)( rand()/(float) RAND_MAX );
}

void init_collisions(bool *collisions, uint max)
{
    for (uint i=0; i < max; ++i)
        collisions[i] = 0;
}


int main(int argc, char** argv)
{
    uint time_stamp = 0, time_max;
    float pressure = 0;
    int c;
    uint nr_particles = INIT_NO_PARTICLES;
    int horiz_size = BOX_HORIZ_SIZE;
    int vert_size = BOX_VERT_SIZE;
    int my_id;  /* My MPI ID */
    int np;  /* Nr. of processors */
    int total_number_of_particles;
    int dims[2];
    int my_coords[2];
    int right_nbr[2];
    int right_nbr_id;

    // Parse arguments
    while ((c = getopt(argc, argv, "hn:x:y:v")) != -1)
    {
        switch (c)
        {
            case 'h':
                help(argv[0]);
                exit(0);
            case 'n':
                nr_particles = atoi(optarg);
                if (nr_particles > MAX_NO_PARTICLES)
                {
                    log_error("The provided number exceeds the maximum (%d).",
                            MAX_NO_PARTICLES);
                    help(argv[0]);
                    exit(1);
                }
                break;
            case 'x':
                horiz_size = atoi(optarg);
                break;
            case 'y':
                vert_size = atoi(optarg);
                break;
            case 'v':
                log_enable(true);
                break;
            default:
                help(argv[0]);
                exit(1);
        }
    }
    if (optind == argc)
    {
        log_error("Error: simulation time is missing");
        help(argv[0]);
        exit(1);
    }
    time_max = atoi(argv[optind]);

    log_info("Nr. particles: %u, Box hz. size: %d, Box vt. size: %d, Simulation time: %d.",
            nr_particles, horiz_size, vert_size, time_max);

    /* Init MPI Environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    log_info("#processors:\t%d, my ID:\t%d", np, my_id);

    /**********************************************************************/
    /* The key now is to introduce new topologies and rearrange
     * processors in a 2D matrix */
    compute_grid_dimensions(np, dims);

    MPI_Dims_create( np, 2, dims);
    int periods[2];
    periods[0] = 1;  /* Periodic in rows */
    periods[1] = 0;  /* Non-periodic in columns */

    int reorder = 1; /* Reorder allowed */

    /* Create the grid topology */
    MPI_Comm grid_comm;

    /* And create the Cart out of the grid object */
    MPI_Cart_create( MPI_COMM_WORLD, 2, dims, periods, reorder, &grid_comm );

    MPI_Cart_get( grid_comm, 2, dims, periods, my_coords );

    MPI_Cart_rank( grid_comm, my_coords, &my_id );

    right_nbr[0] = my_coords[0] + 1;
    right_nbr[1] = my_coords[1];

    MPI_Cart_rank ( grid_comm, right_nbr, & right_nbr_id );
    /**********************************************************************/

    total_number_of_particles = nr_particles * np;
    if (my_id == ROOT)
    {
        log_info("Total amount of existent particles in the box:\t%d", total_number_of_particles);
        log_info("(Grid dimensions are: %d x %d)", dims[0], dims[1]);
    }
    log_info("[ID=%d] My process rank cartesian coordinates: (x=%d, y=%d).", my_id, my_coords[0], my_coords[1]);

    /* Declare our MPI datatypes */
    DECLARE_MPI_COORDINATE_DATATYPE(cord_t, cord_mpi_t);
    DECLARE_MPI_PARTICLE_COORDINATE_DATATYPE(pcord_t, pcord_mpi_t);
    DECLARE_MPI_PARTICLE_DATATYPE(particle_t, particle_mpi_t, pcord_mpi_t);

    /* Initialize */
    // 1. set the walls
    cord_t wall;
    wall.y0 = wall.x0 = 0;
    wall.x1 = horiz_size;
    wall.y1 = vert_size;

    /* Initialize the random seed */
    srand( time(NULL) + 1234 );

    // 2. allocate particle bufer and initialize the particles
    // Note: every processor will have initialized nr_particles, so the total
    // number of particles in the box shall be nr_particles * np
    pcord_t * particles = (pcord_t*) malloc (nr_particles * sizeof(pcord_t) );
    bool * collisions = (bool *) malloc (nr_particles * sizeof(bool) );

    float r, a;
    for (uint i = 0; i < nr_particles; i++)
    {
        // initialize random position
        particles[i].x = wall.x0 + rand1() * horiz_size;
        particles[i].y = wall.y0 + rand1() * vert_size;

        // initialize random velocity
        r = rand1() * MAX_INITIAL_VELOCITY;
        a = rand1() * 2 * PI;
        particles[i].vx = r * cos(a);
        particles[i].vy = r * sin(a);
    }

    uint p, pp;

    /* Start simulation */
    double starttime, endtime;
    starttime = MPI_Wtime();

    /* Main loop */
    for (time_stamp = 0; time_stamp < time_max; time_stamp++) // for each time stamp
    {
        init_collisions(collisions, nr_particles);

        for (p = 0; p < nr_particles; p++) // for all particles
        {
            if (collisions[p]) continue;

            /* check for collisions */
            for (pp = p + 1; pp < nr_particles; pp++)
            {
                if (collisions[pp]) continue;

                float t = collide(&particles[p], &particles[pp]);
                if (t != -1) // collision
                {
                    collisions[p] = collisions[pp] = 1;
                    interact(&particles[p], &particles[pp], t);
                    break; // only check collision of two particles
                }
            }
        }

        // move particles that has not collided with another
        for (p = 0; p < nr_particles; p++)
        {
            if (!collisions[p])
            {
                feuler(&particles[p], 1);

                /* check for wall interaction and add the momentum */
                pressure += wall_collide(&particles[p], wall);
            }
        }
    }

    endtime = MPI_Wtime();

    /* Synchronization point */
    MPI_Barrier(MPI_COMM_WORLD);

    if (my_id == ROOT)
    {
        /* Don't use log.h functions: we want to always output something */
        printf("Average pressure = %f, elapsed time = %.2f secs\n",
                pressure / (WALL_LENGTH * time_max),
                endtime - starttime);
    }

    free(particles);
    free(collisions);

    /* Finalyze MPI running environment */
    MPI_Finalize();

    return 0;
}

