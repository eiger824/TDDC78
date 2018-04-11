#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#include "ppmio.h"
#include "thresfilter.h"

int main (int argc, char ** argv) {
    int xsize, ysize, colmax;
    pixel src[MAX_PIXELS];
    struct timespec stime, etime;

    /* Take care of the arguments */

    if (argc != 3) {
        fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
        exit(1);
    }

#ifdef WITH_MPI
    int p;
    int my_id;
    int elems_per_node;
    uint avg;
    pixel * myarr;
    pixel * myarr2;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    printf("Nprocs:\t%d, my_id:\t%d\n", p, my_id);
#endif

    /* read file */
    if(read_ppm (argv[1], &xsize, &ysize, &colmax, (char *) src) != 0)
        exit(1);

    if (colmax > 255) {
        fprintf(stderr, "Too large maximum color-component value\n");
        exit(1);
    }

    int max_size = xsize * ysize;

    if (my_id == 0)
        printf("Has read the image (%d pixels), calling filter\n", max_size);

#ifdef WITH_MPI
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

    /* Root node */
    if (my_id == 0)
    {
        /* TODO: find out if SIZE is not divisible by p */
        elems_per_node = max_size / p;
        printf("[ROOT] Every processing node will process %d elements.\n" ,elems_per_node);
    }

#endif

#ifdef WITH_MPI
    double starttime, endtime;
    starttime = MPI_Wtime();
#endif

#ifdef WITH_OPENMP
    clock_gettime(CLOCK_REALTIME, &stime);
#endif

#ifdef WITH_MPI

    /* Allocate a processing array for every node */
    myarr = (pixel *) malloc (elems_per_node * sizeof(pixel));

    /* Call MPI_Scatter to compute the sum of all pixels */
    MPI_Scatter((void *)src, max_size, pixel_t_mpi, (void *)myarr, elems_per_node, pixel_t_mpi, 0, MPI_COMM_WORLD);
    uint mysum = get_px_sum(myarr, elems_per_node);

    MPI_Reduce( &mysum, &avg, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );

    /* Now `avg` has the pixel average that we can use in the threshold function */
    if (my_id == 0)
        avg /= max_size;


    myarr2 = (pixel *) malloc (elems_per_node * sizeof(pixel));

    MPI_Scatter( (void *)src, max_size, pixel_t_mpi, (void *) myarr2, elems_per_node, pixel_t_mpi, 0, MPI_COMM_WORLD);
#endif

    /* Call the filtering function */
    thresfilter(myarr2, elems_per_node, avg);

#ifdef WITH_MPI
    /* Gather the results */
    MPI_Gather(myarr2, elems_per_node, pixel_t_mpi, src, max_size, pixel_t_mpi, 0, MPI_COMM_WORLD);
#endif

#ifdef WITH_MPI
    endtime = MPI_Wtime();
    printf("Filtering took: %f secs\n", endtime - starttime);
#endif

#ifdef WITH_OPENMP
    clock_gettime(CLOCK_REALTIME, &etime);
    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
            1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;
#endif

    if (my_id == 0)
    {
        /* write result */
        printf("Writing output file\n");

        if(write_ppm (argv[2], xsize, ysize, (char *)src) != 0)
            exit(1);
    }

#ifdef WITH_MPI
    /* Free allocated buffers */
    free(myarr);
    free(myarr2);

    /* Finalize MPI running environment */
    MPI_Finalize();
#endif

    return(0);
}
