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
#include "blurfilter.h"
#include "gaussw.h"

#define MAX_RAD 1000

#ifdef WITH_PTHREADS

/* Mutex object to use for critical section locks */
static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

/* Thread data struct */
typedef struct thread_data_struct
{
    uint thread_id;       /* Thread-ID allocated */
    int x_start;          /* Where at X coordinate to start */
    int y_start;          /* Where at Y coordinate to start */
    int x_end;            /* Where at X coordinate to end */
    int y_end;            /* Where at Y coordinate to end */
    int x_size;           /* How many horizontal pixels */
    int radius;           /* The radius to use */
    double * w;           /* The gaussian weight array */
    pixel * data_to_process;/* Pointer containing the starting address of
                              the image */
} tdata_t;

/* Threads call this starting routine to process rows */
void * blurf_x_wrapper(void * data)
{
    tdata_t * thread_data = (tdata_t *) data;
    uint id = thread_data->thread_id;
    int xstart = thread_data->x_start;
    int ystart = thread_data->y_start;
    int xend = thread_data->x_end;
    int yend = thread_data->y_end;
    int xsize = thread_data->x_size;
    int radius = thread_data->radius;
    double * weights = thread_data->w;
    pixel * image = thread_data->data_to_process;

    pthread_mutex_lock(&mutex);
    printf("[T%u] hz (xs=%d,ys=%d,xe=%d,ye=%d])\n", id,
            xstart,ystart,xend,yend);
    pthread_mutex_unlock(&mutex);

    // Call the filter in the horizontal direction 
    blurfilter_x(xstart, ystart, xend, yend, xsize, image, radius, weights);

    pthread_exit(NULL);
    return NULL;
}

/* Threads call this starting routine to process columns */
void * blurf_y_wrapper(void * data)
{
    tdata_t * thread_data = (tdata_t *) data;
    uint id = thread_data->thread_id;
    int xstart = thread_data->x_start;
    int ystart = thread_data->y_start;
    int xend = thread_data->x_end;
    int yend = thread_data->y_end;
    int xsize = thread_data->x_size;
    int radius = thread_data->radius;
    double * weights = thread_data->w;
    pixel * image = thread_data->data_to_process;

    pthread_mutex_lock(&mutex);
    printf("[T%u] vt (xs=%d,ys=%d,xe=%d,ye=%d])\n", id,
            xstart,ystart,xend,yend);
    pthread_mutex_unlock(&mutex);

    // Call the filter in the horizontal direction 
    blurfilter_y(xstart, ystart, xend, yend, xsize, image, radius, weights);

    pthread_exit(NULL);
    return NULL;
}

#endif

int main (int argc, char ** argv) {
    int radius;
    int xsize, ysize, colmax;
    pixel src[MAX_PIXELS];
    double w[MAX_RAD];
#ifdef WITH_PTHREADS
    uint nr_threads;
    struct timespec stime, etime;
#endif

    /* Take care of the arguments */

#ifdef WITH_MPI
    if (argc != 4)
    {
        fprintf(stderr, "Usage: %s radius infile outfile\n", argv[0]);
        exit(1);
    }
#endif
#ifdef WITH_PTHREADS
    if (argc != 5)
    {
        fprintf(stderr, "Usage: %s radius infile outfile nr_threads\n", argv[0]);
        exit(1);
    }
    nr_threads = atoi(argv[4]);
#endif

    radius = atoi(argv[1]);

    if (radius > MAX_RAD || radius < 1)
    {
        fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
        exit(1);
    }

    int my_id;
#ifdef WITH_MPI
    int p;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    printf("Nprocs:\t%d, my_id:\t%d\n", p, my_id);
#endif

    /* read file */
    if (read_ppm (argv[2], &xsize, &ysize, &colmax, (char *) src) != 0)
        exit(1);

    if (colmax > 255)
    {
        fprintf(stderr, "Too large maximum color-component value\n");
        exit(1);
    }

    printf("Has read the image, generating coefficients\n");

    /* filter */
    get_gauss_weights(radius, w);
    uint i;
    for (i = 0; i < radius; ++i)
        printf("%.2f ", w[i]);

    printf("\nCalling filter\n");

    /* Different ways of measuring time, according to which technique to use */
#ifdef WITH_MPI
    double starttime, endtime;
    starttime = MPI_Wtime();
#endif

#ifdef WITH_PTHREADS
    pthread_t threads[nr_threads];
    /* Calculate how many rows/colums each thread needs to process */
    /* The first part of the filter runs "horizontally", i.e., 
     * the gaussian weights are applied to horizontally contiguous
     * pixels. Hence, by taking `ysize` and dividing it by `nr_threads`
     * we will obtain how many horizontal blocks each thread will take
     * care of simultaneously. The same applies for the "vertical"
     * filter, where each thread will process a group of columns specified
     * by `nr_threads` and `size`*/
    int hz_block_count = ysize / nr_threads;
    int vt_block_count = xsize / nr_threads;

    tdata_t t[nr_threads];
    clock_gettime(CLOCK_REALTIME, &stime);
    /* Now, create the threads to work on the horizontal filter */
    for (my_id = 1; my_id <= nr_threads; ++my_id)
    {
        t[my_id-1].thread_id = my_id;
        t[my_id-1].x_start = 0;
        t[my_id-1].y_start = (my_id-1)*hz_block_count;
        t[my_id-1].x_end = xsize;
        t[my_id-1].y_end = my_id*hz_block_count; // We want divide the work horizontally
        t[my_id-1].x_size = xsize;
        t[my_id-1].radius = radius;
        t[my_id-1].w = w;
        t[my_id-1].data_to_process = src;

        if (pthread_create(&threads[my_id-1], NULL, blurf_x_wrapper, (void *)&t[my_id-1]) != 0)
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
    printf("[T0] Applied horizontal filter, now applying vertical\n");
    /* Finally, perform the same process for the "vertical" filter */
    for (my_id = 1; my_id <= nr_threads; ++my_id)
    {
        t[my_id-1].thread_id = my_id;
        t[my_id-1].x_start = (my_id-1)*vt_block_count;
        t[my_id-1].y_start = 0; 
        t[my_id-1].x_end = my_id*vt_block_count;
        t[my_id-1].y_end = ysize; // We want divide the work horizontally
        t[my_id-1].x_size = xsize;
        t[my_id-1].radius = radius;
        t[my_id-1].w = w;
        t[my_id-1].data_to_process = src;

        if (pthread_create(&threads[my_id-1], NULL, blurf_y_wrapper, (void *)&t[my_id-1]) != 0)
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
#endif

#ifdef WITH_MPI
    /* Apply the filter */
    blurfilter_x(xsize, ysize, src, radius, w);
    blurfilter_y(xsize, ysize, src, radius, w);

    endtime = MPI_Wtime();
    printf("Filtering took: %f secs\n", endtime - starttime);
#endif

#ifdef WITH_PTHREADS
    clock_gettime(CLOCK_REALTIME, &etime);
    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
            1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;
#endif


    /* write result */
    printf("Writing output file\n");

    if(write_ppm (argv[3], xsize, ysize, (char *)src) != 0)
        exit(1);

    /* Finalize MPI running environment */
#ifdef WITH_MPI
    MPI_Finalize();
#endif

    return 0;
}
