#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <pthread.h>

#include "ppmio.h"
#include "thresfilter.h"

#define     ROOT    0

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

int main (int argc, char ** argv) {
    int xsize, ysize, colmax;
    pixel src[MAX_PIXELS];

    /* Take care of the arguments */
    if (argc != 4)
    {
        fprintf(stderr, "Usage: %s infile outfile nr_threads\n", argv[0]);
        exit(1);
    }

    int elems_per_node = 0;
    int max_size;
    uint avg;
    int my_id;
    /* read file */
    if(read_ppm (argv[1], &xsize, &ysize, &colmax, (char *) src) != 0)
        exit(1);

    if (colmax > 255) {
        fprintf(stderr, "Too large maximum color-component value\n");
        exit(1);
    }
    max_size = xsize * ysize;
    printf("[ROOT] Has read the image (%d pixels), calling filter\n", max_size);

    struct timespec stime, etime;

    int nr_threads = atoi(argv[3]);
    elems_per_node = max_size / nr_threads;

    /* Initalize g_sum_partials array */
    g_sum_partial = (uint *) malloc(nr_threads * sizeof(uint));

    /* Pthread array to use */
    pthread_t threads[nr_threads];
    tdata_t t[nr_threads];
    /* First, all compute the partial sums */
    /* Start by 1, we want the main thread to be 0 - like 'root' in MPI */
    clock_gettime(CLOCK_REALTIME, &stime);
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
    clock_gettime(CLOCK_REALTIME, &etime);
    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
            1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;
    // Free global array of partial sums
    free(g_sum_partial);

    /* write result */
    printf("[ROOT] Writing output file\n");

    if (write_ppm (argv[2], xsize, ysize, (char *)src) != 0)
        exit(1);

    return 0;
}
