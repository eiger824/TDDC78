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
    printf("-n nr\tSet number of particles\n");
    printf("-x size\tSet box horizontal size\n");
    printf("-y size\tSet box vertical size\n");
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

    /* Declare our MPI datatypes */
    DECLARE_MPI_COORDINATE_DATATYPE(cord_t, cord_mpi_t)
    DECLARE_MPI_PARTICLE_COORDINATE_DATATYPE(pcord_t, pcord_mpi_t)
    DECLARE_MPI_PARTICLE_DATATYPE(particle_t, particle_mpi_t, pcord_mpi_t)

    /* Initialize */
    // 1. set the walls
    cord_t wall;
    wall.y0 = wall.x0 = 0;
    wall.x1 = horiz_size;
    wall.y1 = vert_size;


    // 2. allocate particle bufer and initialize the particles
    pcord_t * particles = (pcord_t*) malloc(nr_particles * sizeof(pcord_t) );
    bool * collisions = (bool *) malloc (nr_particles * sizeof(bool) );

    srand( time(NULL) + 1234 );

    float r, a;
    for (uint i = 0; i < nr_particles; i++)
    {
        // initialize random position
        particles[i].x = wall.x0 + rand1()*horiz_size;
        particles[i].y = wall.y0 + rand1()*vert_size;

        // initialize random velocity
        r = rand1()*MAX_INITIAL_VELOCITY;
        a = rand1()*2*PI;
        particles[i].vx = r*cos(a);
        particles[i].vy = r*sin(a);
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

    /* Don't use log.h functions: we want to always output something */
    printf("Average pressure = %f, elapsed time = %.2f secs\n",
            pressure / (WALL_LENGTH * time_max),
            endtime - starttime);

    free(particles);
    free(collisions);

    /* Finalyze MPI running environment */
    MPI_Finalize();

    return 0;
}

