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
#include "utils.h"


//Feel free to change this program to facilitate parallelization.

void help(const char * program)
{
    printf("Usage: %s [h|n|x|y] <sim-time>\n", program);
    printf("-h\tPrint this help and exit\n");
    printf("-n nr\tSet number of particles each processor will handle\n\t(default: %d, max: %d)\n",
            INIT_NO_PARTICLES, MAX_NO_PARTICLES);
    printf("-N nr\tSet the total number of particles in the box\n");
    printf("\t(note: this number should not exceed %d * <nr-available-cores>)\n", MAX_NO_PARTICLES);
    printf("-p\tPrint a scaled representation of the final box\n");
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
    int show_box = 0;
    uint nr_particles = INIT_NO_PARTICLES;
    int horiz_size = BOX_HORIZ_SIZE;
    int vert_size = BOX_VERT_SIZE;
    int my_id;  /* My MPI ID */
    int np;  /* Nr. of processors */
    int total_number_of_particles = -1;
    int dims[2];
    int my_coords[2];
    int horizontal_nbr[2];
    int vertical_nbr[2];
    int horizontal_nbr_id;
    int vertical_nbr_id;
    int mutual_excl_flag = 0;

    // Parse arguments
    while ((c = getopt(argc, argv, "hn:N:px:y:v")) != -1)
    {
        switch (c)
        {
            case 'h':
                help(argv[0]);
                exit(0);
            case 'n':
                if (mutual_excl_flag)
                {
                    log_error("Error: either set the particles per core (-n) or the total amount of particles (-N), not both.");
                    help(argv[0]);
                    exit(1);
                }
                nr_particles = atoi(optarg);
                if (nr_particles > MAX_NO_PARTICLES || nr_particles < INIT_NO_PARTICLES)
                {
                    log_error("The provided number %s (%d).",
                            (nr_particles > MAX_NO_PARTICLES ) ? "exceeds the maximum" : "is below the minimum",
                            (nr_particles > MAX_NO_PARTICLES) ? MAX_NO_PARTICLES : INIT_NO_PARTICLES);
                    help(argv[0]);
                    exit(1);
                }
                mutual_excl_flag = 1;
                break;
            case 'N':
                if (mutual_excl_flag)
                {
                    log_error("Error: either set the particles per core (-n) or the total amount of particles (-N), not both.");
                    help(argv[0]);
                    exit(1);
                }
                total_number_of_particles = atoi(optarg);
                mutual_excl_flag = 1;
                break;
            case 'p':
                show_box = 1;
                break;
            case 'x':
                horiz_size = atoi(optarg);
                break;
            case 'y':
                vert_size = atoi(optarg);
                break;
            case 'v':
                log_enable(INFO);
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

    /* Init MPI Environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    log_info("#processors:\t%d, my ID:\t%d", np, my_id);

    /* Depending on whether the -n or -N flag was provided, the other parameter
     * is calculated out of the provided one*/
    if (total_number_of_particles == -1)
    {
        total_number_of_particles = nr_particles * np;
    }
    else
    {
        nr_particles = total_number_of_particles / np;
        /* We have a requirement of a minimum of 500 particles. Check */
        if (nr_particles < INIT_NO_PARTICLES)
        {
            if (my_id == ROOT)
            {
                log_error("Error: there must be a minimum of %d particles per core.", INIT_NO_PARTICLES);
                help(argv[0]);
            }
            MPI_Finalize();
            exit(1);
        }
    }

    log_info("Nr. particles: %u (%u per core, %u cores), Box hz. size: %d, Box vt. size: %d, Simulation time: %d.",
            total_number_of_particles, nr_particles, np, horiz_size, vert_size, time_max);

    /**********************************************************************/
    /* The key now is to introduce new topologies and rearrange
     * processors in a 2D matrix */
    compute_grid_dimensions(np, dims);

    MPI_Dims_create( np, 2, dims);
    int periods[2];
    periods[0] = 1;  /* Periodic in rows */
    periods[1] = 1;  /* Periodic in columns */

    int reorder = 1; /* Reorder allowed */

    /* Create the grid topology */
    MPI_Comm grid_comm;

    /* And create the Cart out of the grid object */
    MPI_Cart_create( MPI_COMM_WORLD, 2, dims, periods, reorder, &grid_comm );

    MPI_Cart_get( grid_comm, 2, dims, periods, my_coords );

    MPI_Cart_rank( grid_comm, my_coords, &my_id );

    /* Determine my neighbors */
    horizontal_nbr[0] = my_coords[0] + 1;
    horizontal_nbr[1] = my_coords[1];
    vertical_nbr[0] = my_coords[0];
    vertical_nbr[1] = my_coords[1] + 1;

    MPI_Cart_rank ( grid_comm, horizontal_nbr, & horizontal_nbr_id );
    MPI_Cart_rank ( grid_comm, vertical_nbr, & vertical_nbr_id );
    /**********************************************************************/

    if (my_id == ROOT)
    {
        log_info("Total amount of existent particles in the box:\t%d", total_number_of_particles);
        log_info("(Grid dimensions are: %d x %d)", dims[0], dims[1]);
    }
    log_info("[ID=%d] (x=%d, y=%d). My horizontal neighbor [ID=%d] (x=%d ,y=%d). My vertical neighbor [ID=%d] (x=%d, y=%d)."
            , my_id, my_coords[0], my_coords[1]
            , horizontal_nbr_id, horizontal_nbr[0] , horizontal_nbr[1]
            , vertical_nbr_id, vertical_nbr[0] , vertical_nbr[1]);

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

    // 2. allocate particle buffers and initialize the particles
    pcord_t * particles   = (pcord_t *) malloc (nr_particles * sizeof(pcord_t) );
    bool * collisions  = (bool *)    malloc (nr_particles * sizeof(bool) );
    
    /* The particles will now be generated in each processor's grid region */
    float r, a;
    uint xstart, xend;
    uint ystart, yend;
    cord_t * limits;
   
    limits = get_my_grid_boundaries(horiz_size, vert_size, my_coords, dims);
    xstart  = limits->x0; 
    xend    = limits->x1;
    ystart  = limits->y0; 
    yend    = limits->y1;

    print_limits(my_id, limits);

    for (uint i = 0; i < nr_particles; i++)
    {
        // initialize random position in my assigned grid
        particles[i].x = xstart + rand1() * xend;
        particles[i].y = ystart + rand1() * yend;

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

    /* IDEA: Every core will be assigned a grid region, determined by
     * `my_coords`. Whenever a certain particle "hits" the boundary
     * between grid regions, i.e., a particle may be "changing" regions,
     * we need to check if my neighbor has another particle that will
     * collide with our particle. Which neighbor will be determined by
     * which boundary our particle is traversing to: either our horizontal
     * neighbor or our vertical neighbour. If a collision is bound to
     * happen at that particular simulation time, we need to interact
     * both particles. */

    /* Main loop */
    for (time_stamp = 0; time_stamp < time_max; time_stamp++) // for each time stamp
    {
        log_debug("At time stamp %d", time_stamp);
        /* Sett the whole collisions array to false */
        init_collisions(collisions, nr_particles);

        /* Check for collisions within the same grid region*/
        for (p = 0; p < nr_particles; p++) // for all particles
        {
            if ( !collisions[p] )
            {
                for (pp = p + 1; pp < nr_particles; pp++)
                {
                    if ( !collisions[pp] )
                    {
                        float t = collide(&particles[p], &particles[pp]);
                        if (t != -1) // collision
                        {
                            collisions[p] = collisions[pp] = 1;
                            interact(&particles[p], &particles[pp], t);

                            log_debug("[ID=%d] Particles p%u(%.2f,%.2f) and p%u(%.2f,%.2f) will collide in t=%.2f",
                                    my_id, p, particles[p].x, particles[p].y,
                                    pp, particles[pp].x, particles[pp].y, t);
                            break; // only check collision of two particles
                        }
                    }
                }
            }
        }

        /* Move particles that has not collided with another */
        for (p = 0; p < nr_particles; p++)
        {
            if ( !collisions[p] )
            {
                /* Move the particle */
                feuler(&particles[p], 1);

                /* Now we want to check if a particle in my grid region is
                 * close to colliding with another particle in my neighbor's
                 * region */
                if (is_particle_in_grid_boundary(&particles[p], limits)) 
                {
                    /* Then I want to ask my neighbor(s) if they have a similar
                     * particle in our shared boundary in order to interact them.
                     *
                     * Possible scenario:
                     *
                     * MPI_Bcast( ... send the particle here ... );
                     *
                     * Other nodes receive this particle.
                     * They all check if they have a particle in the neighboring
                     * limit.
                     * */
                }

                /* check for wall interaction and add the momentum */
                pressure += wall_collide(&particles[p], wall);
            }
        }
    }

    endtime = MPI_Wtime();

    /* Synchronization point */
    MPI_Barrier(MPI_COMM_WORLD);

    /* Gather all particles if print (-p) option enabled, not mandatory */
    pcord_t * all; 
    if (show_box)
    {
        all = (pcord_t * ) malloc(sizeof(pcord_t) * total_number_of_particles);
        MPI_Gather(particles, nr_particles, pcord_mpi_t,
                all, nr_particles, pcord_mpi_t,
                ROOT, MPI_COMM_WORLD);
    }
    if (my_id == ROOT)
    {
        /* Don't use log.h functions: we want to always output something */
        printf("Average pressure = %f, elapsed time = %.6f secs\n",
                pressure / (WALL_LENGTH * time_max),
                endtime - starttime);
        if (show_box)
        {
            print_box(horiz_size, vert_size, all, total_number_of_particles, dims);
        }
    }
    free(particles);
    free(collisions);

    /* Finalyze MPI running environment */
    MPI_Finalize();

    return 0;
}

