#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#ifdef WITH_MPI
#include <mpi.h>
#else
#include <pthread.h>
#endif

#include "ppmio.h"
#include "thresfilter.h"

#define     ROOT    0

#ifdef WITH_PTHREADS

/* Global sum of pixels */
static uint g_sum = 0;
/* Partial sum of pixels */
static uint * g_sum_partial = NULL;
/* Mutex object to use for critical section locks */
static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

/* Thread data struct  */
typedef struct thread_data_struct 
{
    uint thread_id;  /* Thread-ID allocated */
    uint nr_elems;   /* How many elems will I process? */
    uint px_avg;     /* The pixel intensity average */
    void * data_to_process;  /* Pointer containing the address of
                                the segment of the array to process */
} tdata_t;

void * get_px_sum_wrapper(void * data)
{
    tdata_t * thread_data = (tdata_t * ) data;
    uint id = thread_data->thread_id;
    uint elems = thread_data->nr_elems;
    pixel * array = (pixel *) thread_data->data_to_process;

    uint sum = get_px_sum(array, elems);
    g_sum_partial[id-1] = sum;

    pthread_mutex_lock(&mutex);
    printf("[T%u] Calculating my partial sum:\t%u\n", id, sum);
    pthread_mutex_unlock(&mutex);

    pthread_exit(NULL);
    return NULL;
}

void * thresfilter_wrapper(void * data)
{
    tdata_t * thread_data = (tdata_t * ) data;
    uint id = thread_data->thread_id;
    uint elems = thread_data->nr_elems;
    uint avg = thread_data->px_avg;
    pixel * array = (pixel *) thread_data->data_to_process;

    pthread_mutex_lock(&mutex);
    printf("[T%u] Applying filter on my image segment\n", id);
    pthread_mutex_unlock(&mutex);

    thresfilter(array, elems, avg);

    pthread_exit(NULL);
    return NULL;
}

uint get_average(uint nr_threads, uint max_pixels)
{
    uint avg = 0;
    uint i;
    for (i = 0; i < nr_threads; ++i)
    {
        avg += g_sum_partial[i];
    }
    avg /= max_pixels;
    return avg;
}

#endif

int main (int argc, char ** argv) {
    int xsize, ysize, colmax;
    pixel src[MAX_PIXELS];
    struct timespec stime, etime;

    /* Take care of the arguments */
#ifdef WITH_MPI
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
        exit(1);
    }
#endif
#ifdef WITH_PTHREADS
    if (argc != 4)
    {
        fprintf(stderr, "Usage: %s infile outfile nr_threads\n", argv[0]);
        exit(1);
    }
#endif

    int elems_per_node = 0;
    int max_size;
    uint avg;
    int my_id;
#ifdef WITH_MPI
    int p;
    pixel * myarr;
    pixel * myarr2;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    printf("Nprocs:\t%d, my_id:\t%d\n", p, my_id);
#endif

#ifdef WITH_MPI
    if (my_id == ROOT)
    {
#endif
        /* read file */
        if(read_ppm (argv[1], &xsize, &ysize, &colmax, (char *) src) != 0)
            exit(1);

        if (colmax > 255) {
            fprintf(stderr, "Too large maximum color-component value\n");
            exit(1);
        }
        max_size = xsize * ysize;
        printf("[ROOT] Has read the image (%d pixels), calling filter\n", max_size);
        /* Once root has read and computed all pixels, send the max pixels to all processes */
#ifdef WITH_MPI
        unsigned i;
        for (i = 1; i < p; ++i)
            MPI_Send(&max_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
    else
        MPI_Recv(&max_size, 1, MPI_INT, ROOT, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

    /************* First create our own MPI datatype *************/
    pixel item;
    MPI_Datatype pixel_t_mpi; // MPI type to commit

    int block_lengths [] = {1, 1, 1};
    MPI_Datatype block_types[] = { MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR };
    MPI_Aint start, displ[3];

    MPI_Address( &item, &start );
    MPI_Address( &item.r, &displ[0]);
    MPI_Address( &item.g, &displ[1]);
    MPI_Address( &item.b, &displ[2]);

    displ[0] -= start;
    displ[1] -= start;
    displ[2] -= start;

    MPI_Type_struct( 3, block_lengths, displ, block_types, &pixel_t_mpi );
    /* Commit the newly defined type */
    MPI_Type_commit ( &pixel_t_mpi );

    /************ END OF MPI DATATYPE CREATION ************/

    elems_per_node = max_size / p;

    /* Root node */
    if (my_id == ROOT)
    {
        /* TODO: find out if SIZE is not divisible by p */
        printf("[ROOT] Every processing node will process %d elements.\n" ,elems_per_node);
    }

    double starttime, endtime;
    starttime = MPI_Wtime();
#endif

#ifdef WITH_PTHREADS
    clock_gettime(CLOCK_REALTIME, &stime);

    /* PTHread array to use */
    int nr_threads = atoi(argv[3]);
    pthread_t threads[nr_threads];
    elems_per_node = max_size / nr_threads;

    /* Initalize g_sum_partials array */
    g_sum_partial = (uint *) malloc(nr_threads * sizeof(uint));
    tdata_t t[nr_threads];
    /* First, all compute the partial sums */
    /* Start by 1, we want the main thread to be 0 - like 'root' in MPI */
    for (my_id = 1; my_id <= nr_threads; ++my_id)
    {
        t[my_id-1].thread_id = my_id;
        t[my_id-1].nr_elems = elems_per_node;
        t[my_id-1].data_to_process = (void *)(src + (my_id-1)*elems_per_node);
        if (pthread_create(&threads[my_id-1], NULL, get_px_sum_wrapper, (void *)&t[my_id-1]) != 0)
        {
            perror("Error creating thread.");
            exit(1);
        }
    }
    /* Next, wait for them to finish and return */
    for (my_id = 1; my_id <= nr_threads; ++my_id)
    {
        if (pthread_join(threads[my_id-1], NULL) != 0)
        {
            perror("Error joining thread.");
            exit(2);
        }
    }
    /* Now, root (main thread) calculates the average from the partial sums */
    avg = get_average(nr_threads, max_size);
    printf("[ROOT] Average was:\t%u\n", avg);
    /* Finally, generate threads (again) to perform the actual filtering */
    for (my_id = 1; my_id <= nr_threads; ++my_id)
    {
        t[my_id-1].thread_id = my_id;
        t[my_id-1].nr_elems = elems_per_node;
        t[my_id-1].px_avg = avg;
        t[my_id-1].data_to_process = (void *)(src + (my_id-1)*elems_per_node);
        if (pthread_create(&threads[my_id-1], NULL, thresfilter_wrapper, (void *)&t[my_id-1]) != 0)
        {
            perror("Error creating thread.");
            exit(1);
        }
    }
    /* And wait for them to finish and return */
    for (my_id = 1; my_id <= nr_threads; ++my_id)
    {
        if (pthread_join(threads[my_id-1], NULL) != 0)
        {
            perror("Error joining thread.");
            exit(2);
        }
    }
#endif

#ifdef WITH_MPI

    /* Allocate a processing array for every node */
    myarr = (pixel *) malloc (elems_per_node * sizeof(pixel));

    /* Call MPI_Scatter to compute the sum of all pixels */
    MPI_Scatter((void *)src, elems_per_node, pixel_t_mpi, (void *)myarr, elems_per_node, pixel_t_mpi, ROOT, MPI_COMM_WORLD);
    uint mysum = get_px_sum(myarr, elems_per_node);

    MPI_Reduce( &mysum, &avg, 1, MPI_UNSIGNED, MPI_SUM, ROOT, MPI_COMM_WORLD );

    /* Now `avg` has the pixel average that we can use in the threshold function */
    avg /= max_size;

    /* Now we need to send this average to the processes in this group */
    unsigned i;
    if (my_id == ROOT)
        for (i = 1; i < p; ++i)
            MPI_Send(&avg, 1, MPI_UNSIGNED, i, ROOT, MPI_COMM_WORLD);
    else
        MPI_Recv(&avg, 1, MPI_UNSIGNED, ROOT, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

    myarr2 = (pixel *) malloc (elems_per_node * sizeof(pixel));

    MPI_Scatter( (void *)src, elems_per_node, pixel_t_mpi, (void *) myarr2, elems_per_node, pixel_t_mpi, ROOT, MPI_COMM_WORLD);

    /* Call the filtering function */
    thresfilter(myarr2, elems_per_node, avg);

    /* Gather the results */
    MPI_Gather(myarr2, elems_per_node, pixel_t_mpi, src, elems_per_node, pixel_t_mpi, ROOT, MPI_COMM_WORLD);

    /* After this point, we know all group members have entered the barrier
     * i.e., all have processed their part of the image */
    MPI_Barrier(MPI_COMM_WORLD);

    endtime = MPI_Wtime();
    if (my_id == ROOT)
        printf("Filtering took: %f secs\n", endtime - starttime);
#endif

#ifdef WITH_PTHREADS
    clock_gettime(CLOCK_REALTIME, &etime);
    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
            1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;
    // Free global array of partial sums
    free(g_sum_partial);
#endif

#ifdef WITH_MPI
    if (my_id == ROOT)
    {
#endif
        /* write result */
        printf("[ROOT] Writing output file\n");

        if (write_ppm (argv[2], xsize, ysize, (char *)src) != 0)
            exit(1);
#ifdef WITH_MPI
    }

    /* Free allocated buffers */
    free(myarr);
    free(myarr2);

    /* Finalize MPI running environment */
    MPI_Finalize();
#endif

    return 0;
}
