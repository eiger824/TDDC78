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

/* This array will hold the ranks of my neighbors */
/*                     L    U   R   B  */
static int my_neighbors[4] = {0,   0,  0,  0};

//Feel free to change this program to facilitate parallelization.
void help(const char * program)
{
    printf("Usage: %s [h|n|x|y] <sim-time>\n", program);
    printf("-h\tPrint this help and exit\n");
    printf("-n nr\tSet number of particles each processor will handle\n\t(default: %d, max: %d)\n",
            INIT_NO_PARTICLES, MAX_NO_PARTICLES);
    printf("-N nr\tSet the total number of particles in the box\n");
    printf("\t(note: this number should not exceed %d * <nr-available-cores>)\n", MAX_NO_PARTICLES);
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

static nloc_t get_nbr_dir(int rank)
{
    for (nloc_t nbr = LEFT; nbr <= BOTTOM; ++nbr)
    {
        if (my_neighbors[nbr] == rank)
            return nbr;
    }
    return -1;
}

static void log_if_0(const char * fmt, ...)
{
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if (id == 0)
    {
        va_list args;
        va_start(args, fmt);
        fprintf(stdout, "[ID=0] ");
        vfprintf(stdout, fmt, args);
        va_end(args);
    }
}

static void print_lists(dll_t ** lists, int id)
{
    printf("[ID=%d] My send lists:", id);
    dll_t * left_list = *(lists + LEFT);
    dll_t * top_list = *(lists + TOP);
    dll_t * right_list = *(lists + RIGHT);
    dll_t * bottom_list = *(lists + BOTTOM);
    printf("L(%d), T(%d), R(%d), B(%d)\n",
            left_list->count, top_list->count, right_list->count, bottom_list->count);
}

int main(int argc, char** argv)
{
    uint time_stamp = 0, time_max;
    float pressure = 0;
    float global_pressure = 0;
    int c;
    int verbose = -1;
    uint nr_particles = 1e4;
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
    while ((c = getopt(argc, argv, "hn:N:t:x:y:v")) != -1)
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
                // Jump straight to WARNING if a single '-v' is provided
                if (verbose == -1)
                    verbose = 1;
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
    periods[0] = 0;  /* Non-Periodic in rows */
    periods[1] = 0;  /* Non-Periodic in columns */

    int reorder = 0; /* Reorder not allowed */

    /* Create the grid topology */
    MPI_Comm grid_comm;

    /* And create the Cart out of the grid object */
    MPI_Cart_create( MPI_COMM_WORLD, 2, dims, periods, reorder, &grid_comm );
    MPI_Cart_get( grid_comm, 2, dims, periods, my_coords );

    /* Compute my neigbors */
    int src;
    MPI_Cart_shift(grid_comm, 1, -1, &src, &my_neighbors[LEFT]);
    MPI_Cart_shift(grid_comm, 0, -1, &src, &my_neighbors[TOP]);
    MPI_Cart_shift(grid_comm, 1,  1, &src, &my_neighbors[RIGHT]);
    MPI_Cart_shift(grid_comm, 0,  1, &src, &my_neighbors[BOTTOM]);

    log_info("I am %d, my neighbors are L=%d, T=%d, R=%d, B=%d",
            my_id, my_neighbors[LEFT], my_neighbors[TOP], my_neighbors[RIGHT], my_neighbors[BOTTOM] );
    /**********************************************************************/

    if (my_id == ROOT)
    {
        log_info("Total amount of existent particles in the box:\t%d", total_number_of_particles);
        log_info("(Grid dimensions are: %d x %d)", dims[0], dims[1]);
    }
    log_debug("[ID=%d] (x=%d, y=%d).", my_id, my_coords[0], my_coords[1]);

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

    /* The main list: the particles I will handle */
    dll_t * my_list = dll_init();

    /* The send lists: particles "abandoning" my region, one for each direction */
    dll_t ** my_send_lists = (dll_t ** )  malloc (sizeof *my_send_lists * 4);
    for (int nbr = LEFT; nbr <= BOTTOM; ++nbr)
        my_send_lists[nbr] = dll_init();

    /* The receive arrays: particles "coming" to my region, one from each direction */
    pcord_t ** my_received_particles = (pcord_t ** ) malloc (sizeof *my_received_particles * 4);
    // Don't allocate individual sizes just yet, do it dynamically

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
    pcord_t particle;
    for (uint i = 0; i < nr_particles; i++)
    {
        // initialize random position in my assigned grid
        particle.x = xstart + rand() % (xend - xstart + 1);
        particle.y = ystart + rand() % (yend - ystart + 1);

        // initialize random velocity
        r = rand1() * MAX_INITIAL_VELOCITY;
        a = rand1() * 2 * PI;
        particle.vx = r * cos(a);
        particle.vy = r * sin(a);

        // Append it to the list
        dll_append(my_list, particle);
    }

    dll_node_t * current_node,      * next_node;
    pcord_t    current_particle,  next_particle, particle_to_send;
    int i,j;

    MPI_Status status;

    /* Start simulation */
    double starttime, endtime;
    starttime = MPI_Wtime();
    double a_time, b_time, cs_time;
    float t;
    bool collided = false;
    long count = 0;

	/* Main loop */
    for (time_stamp = 0; time_stamp < time_max; time_stamp++) // for each time stamp
    {
        if (my_id == ROOT)
            log_warning("[ID=%d] At time stamp %d", my_id, time_stamp);

        /* Check for collisions */
        current_node = my_list->head->next;
        while (current_node != my_list->tail)
        {
            current_particle    = current_node->p;

            collided = false;
            next_node = current_node->next;
            while (next_node != my_list->tail)
            {
                next_particle   = next_node->p;

                if ( (t = collide(&current_particle, &next_particle)) != -1)  // Collision
                {
                    /* Interact particles to collide */
                    interact(&current_particle, &next_particle, t);
                    /* Swap collided particle with next on the list since it's a much
                     * faster operation than keep looping on a false condition.
                     * Of course, do it if the particles are valid */
                    dll_swap_nodes(my_list, current_node->next, next_node);
                    current_node = current_node->next;
                    collided = true;
                    break;
                }
                next_node = next_node->next;
            }

            /* Move particle that has not collided with another */
            if (!collided)
            {
                feuler(&current_particle, time_step);

                /* check for wall interaction and add the momentum */
                pressure += wall_collide(&current_particle, wall);

            }
            current_node = current_node->next;
        }

        current_node = my_list->head->next;
        while (current_node != my_list->tail)
        {
            current_particle    = current_node->p;

            nloc_t direction;
            /* Moving particles from my big list to the send list */
            if (is_particle_outside_grid_boundary(
                        &current_particle, limits,
                        &direction))
            {
                // Extract it from my original list
                printf("Size:\t%d ------->", my_list->count);
                current_node = dll_extract(my_list, current_node, &particle_to_send);
                printf("%d\n", my_list->count);

                /* Only send away this particle if this is a valid neighbor */
                if (direction != -1)
                {
                    dll_t * selected_outlist = my_send_lists[direction];
                    // Append it to this send list
                    dll_append(selected_outlist, particle_to_send);
                }
            }
            else
            {
                current_node = current_node->next;
            }
        }
        // At this point, some particles may have disappeared from "my_list" and
        // others may have appeared in the correspondent list of "my_send_lists".
        // Now we want to handle the particles in "my_send_lists" so they can be
        // sent to their respective destinations

        a_time = MPI_Wtime();
        /* Sending the particles in all directions */
        int recv_sizes[4] = {0};
        // nbr: who am I sending to ??
        for (int nbr = LEFT; nbr <= BOTTOM; ++nbr)
        {
            int opposite = (nbr + 2) % 4;
            int src = my_neighbors[nbr];
            int from = my_neighbors[opposite];
            // 1. Convert the current send list to an appropriate C-style array
            dll_t *   tmp_send_list = my_send_lists[nbr];
            pcord_t * tmp_send_array = dll_to_array( tmp_send_list );
            // 2. Send it to the opposite direction (L->R, R->L, T->B, B->T)
            MPI_Send(tmp_send_array, tmp_send_list->count, pcord_mpi_t, src, nbr, grid_comm);
            // 3. Probe in order to check whether a message is waiting for me 
            MPI_Probe(from, nbr, grid_comm, &status);
            // 4. Get the number of bytes on the receiver buffer (could be 0 ..)
            MPI_Get_count(&status, pcord_mpi_t, &recv_sizes[nbr]);
            // 5. Depending on the obtained size, allocate a receive buffer
            if (recv_sizes[nbr] > 0)
            {
                my_received_particles[nbr] = (pcord_t * ) malloc (sizeof (pcord_t) * recv_sizes[nbr]);
            }
            else
            {
                my_received_particles[nbr] = NULL;
            }
            // 6. Receive the actual message that is waiting for me
            MPI_Recv( my_received_particles[nbr] , recv_sizes[nbr], pcord_mpi_t, from, nbr, grid_comm, &status);
        }
        b_time = MPI_Wtime();
        cs_time += (b_time - a_time);

        /* Go through my incoming particles and append them to my main list
           After every time stamp, remove all elems from all send lists as
           well */
        for (int nbr = LEFT; nbr <= BOTTOM; ++nbr)
        {
            for (int pos = 0; pos < recv_sizes[nbr]; ++pos) 
            {
                pcord_t new = my_received_particles[nbr][pos];
                log_debug("[ID:%d] Appending new!", my_id);
                dll_append(my_list, new);
            }
            if (recv_sizes[nbr] > 0)
                free(my_received_particles[nbr]);

            dll_empty( my_send_lists[nbr] );
        }
        MPI_Barrier(grid_comm);
    }

    MPI_Allreduce(&pressure, &global_pressure, 1, MPI_FLOAT, MPI_SUM, grid_comm);
    endtime = MPI_Wtime();

    if (my_id == ROOT)
    {
        double avg = cs_time;
        double total_elapsed = endtime - starttime;
        /* Don't use log.h functions: we want to always output something */
        printf("Average pressure = %f, elapsed time = %.6f secs\n",
                pressure / (WALL_LENGTH * time_max),
                total_elapsed);
        printf("Average time spent on the communications:\t%.6fs (%.2f %%)\n",
                avg, 100 * avg / total_elapsed);
    }

    /* Free stuff */
    dll_destroy(my_list);

    for (uint i = LEFT; i <= BOTTOM; ++i)
        dll_destroy(my_send_lists[i]);
    free(my_send_lists);

/*
    for (uint i = LEFT; i <= BOTTOM; ++i)
        free(my_received_particles[i]);
    free(my_received_particles);
*/

    /* Finalyze MPI running environment */
    log_debug("[ID=%d] Done. Au revoir!", my_id);
    MPI_Finalize();

    return 0;
}

