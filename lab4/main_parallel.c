#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
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
#include "dll.h"

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
    printf("-t step\tTime step to move the particles when they don't collide\n\t(default: %.2f)\n", DEFAULT_TSTEP);
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

int main(int argc, char** argv)
{
    uint time_stamp = 0, time_max;
    float pressure = 0;
    float global_pressure = 0;
    int c;
    int verbose = -1;
    int show_box = 0;
    uint nr_particles = INIT_NO_PARTICLES;
    int horiz_size = BOX_HORIZ_SIZE;
    int vert_size = BOX_VERT_SIZE;
    int my_id;  /* My MPI Id */
    int np;  /* Nr. of processors */
    int total_number_of_particles = -1;
    int dims[2];
    int my_coords[2];
    float time_step = DEFAULT_TSTEP;
    int mutual_excl_flag = 0;

    // Parse arguments
    while ((c = getopt(argc, argv, "hn:N:pt:x:y:v")) != -1)
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
            case 't':
                time_step = atoi(optarg);
                break;
            case 'x':
                horiz_size = atoi(optarg);
                break;
            case 'y':
                vert_size = atoi(optarg);
                break;
            case 'v':
                verbose++;
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

    log_enable(verbose);

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

    /**********************************************************************/

    if (my_id == ROOT)
    {
        log_info("Total amount of existent particles in the box:\t%d", total_number_of_particles);
        log_info("(Grid dimensions are: %d x %d)", dims[0], dims[1]);
    }
	log_info("[ID=%d] (x=%d, y=%d).", my_id, my_coords[0], my_coords[1]);

    /* Declare our MPI datatypes */
    DECLARE_MPI_COORDINATE_DATATYPE(cord_t, cord_mpi_t);
    DECLARE_MPI_PARTICLE_COORDINATE_DATATYPE(pcord_t, pcord_mpi_t);
    DECLARE_MPI_PARTICLE_DATATYPE(particle_t, particle_mpi_t, pcord_mpi_t);

    /* Initialize */
    cord_t wall;
    wall.y0 = wall.x0 = 0;
    wall.x1 = horiz_size;
    wall.y1 = vert_size;

    /* Initialize the random seed */
    srand( time(NULL) + 1234 );

    /* Initalize the doubly linked list */
    dll_t * my_list = dll_init();

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

    // warning mode only
    print_limits(my_id, limits);

    // Needed types
    pcord_t * particle;
    for (uint i = 0; i < nr_particles; i++)
    {
        // Allocate memory for the particle
        particle = (pcord_t * ) malloc (sizeof *particle);
        // initialize random position in my assigned grid
        particle->x = xstart + rand1() * xend;
        particle->y = ystart + rand1() * yend;

        // initialize random velocity
        r = rand1() * MAX_INITIAL_VELOCITY;
        a = rand1() * 2 * PI;
        particle->vx = r * cos(a);
        particle->vy = r * sin(a);

        // Append it to the list
        dll_append(my_list, particle);
    }


    dll_node_t * current_node,      * next_node;
    pcord_t    * current_particle,  * next_particle, * particle_to_send;
    int i,j;
    bool collided;
    int neighbor_coordinates[2];




    dll_t * my_send_list = dll_init();
    //
    // We want to hold a "list of lists" with where each list will hold
    // the particles to be sent to each processor


    pcord_t ** travelling_particles = (pcord_t ** ) malloc (sizeof *travelling_particles * np);
    /* Allocate place for at most "nr_particles" particles that leave a specific region */
    for (unsigned i = 0; i < np; ++i)
        *(travelling_particles + i) = (pcord_t * ) malloc (sizeof(pcord_t) * nr_particles);

    /* Start simulation */
    double starttime, endtime;
    starttime = MPI_Wtime();

	/* Main loop */
    for (time_stamp = 0; time_stamp < time_max; time_stamp++) // for each time stamp
    {
        log_debug("At time stamp %d", time_stamp);

        /* Check for collisions */
        for (i = 0; i < my_list->count; ++i)
        {
            current_node        = dll_at(my_list, i);
            current_particle    = current_node->p;

            collided = false;

            for (j = i + 1; j < my_list->count; ++j)
            {
                next_node       = dll_at(my_list, j);
                next_particle   = next_node->p;

                float t = collide(current_particle, next_particle);
                if (t != -1)  // Collision
                {
                    collided = true;
                    interact(current_particle, next_particle, t);
                }
            }

            /* Move particle that has not collided with another */
            if (!collided)
                feuler(current_particle, time_step);

            /* check for wall interaction and add the momentum */
            pressure += wall_collide(current_particle, wall);

            /* Moving particles from my big list to the send list */
            if (is_particle_outside_grid_boundary(
                        current_particle, limits,
                        horiz_size, vert_size,
                        dims, neighbor_coordinates)) 
            {
				// Determine what is the neigbor rank
				int nbr_rank;
				MPI_Cart_rank(grid_comm, neighbor_coordinates, &nbr_rank);

                log_debug("Transfer: %d\t->\t%d", my_id, nbr_rank);

                // Extract it from my original list
                particle_to_send = (pcord_t * ) malloc (sizeof *particle_to_send);
                dll_extract(my_list, current_node, particle_to_send);

                // Append it to my send list
                dll_append(my_send_list , particle_to_send);
            }
        }
        // At this point, some particles may have disappeared from "my_list" and
        // others may have appeared in the correspondent list of "my_send_lists".
        // Now we want to handle the particles in "my_send_lists" so they can be
        // sent to their respective destinations

        /* Convert my send-list to an array of coordinates */
        pcord_t * my_particles = dll_to_array(my_send_list);
        memcpy(*(travelling_particles + my_id), my_particles, sizeof(pcord_t) * my_send_list->count);

        /* At this point, the matrix "travelling particles should already be filled */
        MPI_Bcast(*(travelling_particles + my_id), my_send_list->count, pcord_mpi_t, my_id, MPI_COMM_WORLD);

        /* Synchronization point, needed for the next operation */
        MPI_Barrier(MPI_COMM_WORLD);

        /* Loop through the travelling particle matrix and just append those particles
         * that belong to my grid area */
        for (unsigned i = 0; i < np; ++i)
        {
            for (unsigned j = 0; j < nr_particles; ++j)
            {
                pcord_t * current = *(travelling_particles + i) + j;
                if (is_particle_inside_grid_boundary(current, limits))
                {
                    // Append this particle to "my_list"!
                    dll_append(my_list, current);
                    log_debug("[ID=%d] Appended new particle to my list. New size: %d",
                            my_id, my_list->count);
                }
            }
        }

/*
        unsigned sender, receiver;
        MPI_Status status;
        for (sender = 0; sender < np; ++sender)
        {
            for (receiver = 0; receiver < np; ++receiver)
            {
                if (sender != receiver)
                {
                    if (my_id == sender)
                    {
                        MPI_Send(nr_particles_for_every_node, np, MPI_INT, receiver, 1234, MPI_COMM_WORLD);
                    }
                    if (my_id == receiver)
                    {
                        MPI_Recv(*(counts + sender), np, MPI_INT, sender, 1234, MPI_COMM_WORLD, &status);
                    }
                }
            }
        }
*/
        /* ****************************************************************************
         *
         *
         * Remove this section
         *
         * ***************************************************************************/
/*
        if (my_id == ROOT)
        {
            FILE * fp = fopen("counts.txt", "a+");
            if (fp != NULL)
            {
                fprintf(fp, "At time stamp %d, the matrix looked like this:\n", time_stamp);
                for (unsigned i = 0; i < np; ++i)
                {
                    for (unsigned j = 0; j < np; ++j)
                    {
                        fprintf(fp, "%d\t", counts[i][j]);
                    }
                    fprintf(fp, "\n");
                }
                fprintf(fp, "\n\n");
                fclose(fp);
            }
        }

*/
        /* ****************************************************************************
         *
         *
         *
         * ***************************************************************************/
        // Send the particles
        
/*
		for (sender = 0; sender < np; ++sender) // sender -> sending process
		{
			for (receiver = 0; receiver < np; ++receiver) // receiver -> receiving process
			{
				if (sender != receiver) // Dont send to the same node
				{
					if (my_id == sender)
					{
                        log_debug("ID=%d, sending to %d", sender, receiver);
                        dll_t * tmp = *(my_send_lists + receiver);
                        pcord_t * array_to_send = dll_to_array(tmp);
                        if (array_to_send)
                        {
                            log_debug("(ptr '%p' seems legit.)", (void*)array_to_send);
                        }
                        MPI_Send(array_to_send, counts[sender][receiver], pcord_mpi_t,
                                    receiver, 1234, MPI_COMM_WORLD);
					}
					if (my_id == receiver)
					{
                        log_debug("ID=%d, receiving from %d", receiver, sender);
                        MPI_Recv(receive_from_others, counts[sender][receiver], pcord_mpi_t, sender, 1234, MPI_COMM_WORLD, &status);
                        // At this point, "receive_from_others" contains particles that I should append to my main particle list
                        for (unsigned i = 0; i < counts[sender][receiver]; ++i)
                        {
                            pcord_t * new = (pcord_t * ) malloc(sizeof *new);
                            *new = *(receive_from_others + i);
                            dll_append(my_list, new);
                        }
					}
				}
			}
		}
*/
        
        /* After every time stamp, remove all elems from all send lists */
        dll_empty(my_send_list);
        /* And restore the counts */
//         for (unsigned i = 0; i < np; ++i)
//         {
//             for (unsigned j = 0; j < np; ++j)
//             {
//                 counts[i][j] = 0;
//             }
//         }

    }

    MPI_Allreduce(&pressure, &global_pressure, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    endtime = MPI_Wtime();

    /* Synchronization point */
    MPI_Barrier(MPI_COMM_WORLD);

    if (my_id == ROOT)
    {
        /* Don't use log.h functions: we want to always output something */
        printf("Average pressure = %f, elapsed time = %.6f secs\n",
                pressure / (WALL_LENGTH * time_max),
                endtime - starttime);
    }

    /* Free stuff */
#if 0
    dll_destroy(my_list);
    for (unsigned i = 0; i < np; ++i)
    {
        void * ptr1 = (void * ) *(counts + i);
        void * ptr2 = (void * ) *(my_send_lists + i);
//         free(*(counts + i));
//         dll_destroy(*(my_send_lists + i));
        printf("Freeing %p and %p\n", ptr1, ptr2);
        free(ptr1);
        free(ptr2);
    }
    free(counts);
    free(my_send_lists);
#endif

    /* Finalyze MPI running environment */
    log_debug("[ID=%d] Done. Au revoir!", my_id);
    MPI_Finalize();

    return 0;
}

