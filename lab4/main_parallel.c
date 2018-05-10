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

    dll_node_t * current_node,      * next_node;
    pcord_t    * current_particle,  * next_particle, * particle_to_send;
    int i,j;
    bool collided;
    int neighbor_coordinates[2];
    dll_t * my_send_list = dll_init();
    // Allocate a receive buffer big enough to hold at most extra nr_particles
    pcord_t * receive_from_others = (pcord_t * ) malloc (sizeof *receive_from_others * nr_particles);
	// Allocate an array containing how many particles for each node to send (init with zeros)
	int * nr_particles_for_every_node = (int *) calloc (np, sizeof *nr_particles_for_every_node);
	// Same array on the receive nodes
	int * nr_particles_for_every_node_recv = (int *) calloc (np, sizeof *nr_particles_for_every_node_recv);

	/* Main loop */
    for (time_stamp = 0; time_stamp < time_max; time_stamp++) // for each time stamp
    {
        log_debug("At time stamp %d", time_stamp);

        /* At every time stamp, remove all elems from the send list */
        dll_empty(my_send_list);

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
                particle_to_send = (pcord_t * ) malloc (sizeof *particle_to_send);
                // Extract it from my original list
                dll_extract(my_list, current_node, particle_to_send);
                // Append it to my send list
                dll_append(my_send_list, particle_to_send);
				// Determine what is the neigbor rank
				int nbr_rank;
				MPI_Cart_rank(grid_comm, neighbor_coordinates, &nbr_rank);
				// Add up the particles to send for that neighbor
				nr_particles_for_every_node[nbr_rank]++;
            }
        }
        // At this point, some particles may have disappeared from "my_list" and
        // others may have appeared in "my_send_list". Now we want to handle the
        // particles in "my_send_list" so they can be sent to their respective
        // destinations

		int ** matrix = (int ** ) malloc (sizeof *matrix * np);
		for (unsigned i = 0; i < np; ++i)
			*(matrix + i) = (int * ) calloc (np, sizeof(int));

		// At this point, only I know the count of the elements to send
		memcpy(*(matrix + my_id), nr_particles_for_every_node, sizeof(int) * np);

		// Send my row to everyone
		MPI_Bcast(*(matrix + my_id), np, MPI_INT, my_id, MPI_COMM_WORLD);

		// Let's see if matrix works
		for (unsigned i = 0; i < np; ++i)
		{
			for (unsigned j = 0; j < np; ++i)
			{
				printf("%d\t", matrix[i][j]);
			}
			printf("\n");
		}











		//TODO: Replace my_send_list with my_send_array
/*
		unsigned sender, receiver;
		for (sender = 0; sender < np; ++sender) // sender -> sending process
		{
			for (receiver = 0; receiver < np; ++receiver) // receiver -> receiving process
			{
				if (sender != receiver) // Dont send to the same node
				{
					if (my_id == sender)
					{
// 						MPI_Send(my_send_array, nr_elems_to_send, pcord_mpi_t, receiver, 1234, MPI_COMM_WORLD);
					}
					if (my_id == receiver)
					{
						MPI_Status status;
// 						MPI_Recv(my_recv_array, nr_elems_to_send, pcord_mpi_t, sender, 1234, MPI_COMM_WORLD, &status);
					}
				}
			}
		}

*/
    }

    MPI_Allreduce(&pressure, &global_pressure, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    endtime = MPI_Wtime();

    /* Synchronization point */
    MPI_Barrier(MPI_COMM_WORLD);

    /* Gather all particles if print (-p) option enabled, not mandatory */
#if 0
	pcord_t * indiv = dll_to_array(list);
    pcord_t * all; 
    if (show_box)
    {
        all = (pcord_t * ) malloc(sizeof(pcord_t) * total_number_of_particles);
		MPI_Gatherv(indiv , list->count, pcord_mpi_t,
                all, nr_particles, pcord_mpi_t,
                ROOT, MPI_COMM_WORLD);
    }
#endif
    if (my_id == ROOT)
    {
        /* Don't use log.h functions: we want to always output something */
        printf("Average pressure = %f, elapsed time = %.6f secs\n",
                pressure / (WALL_LENGTH * time_max),
                endtime - starttime);
#if 0
        if (show_box)
        {
            print_box(horiz_size, vert_size, all, total_number_of_particles, dims);
        }
#endif
    }
    /* Free both allocated lists */
    dll_destroy(my_list);
    dll_destroy(my_send_list);

    /* Finalyze MPI running environment */
    MPI_Finalize();

    return 0;
}

