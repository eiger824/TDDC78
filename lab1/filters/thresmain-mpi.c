#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <mpi.h>

#include "ppmio.h"
#include "thresfilter.h"

#define     ROOT    0

int main (int argc, char ** argv) {
    int xsize, ysize, colmax;
    pixel * src = (pixel *) malloc (sizeof*src * MAX_PIXELS);

    /* Take care of the arguments */
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
        exit(1);
    }

    int elems_per_node = 0;
    int max_size;
    uint avg;
    int my_id;
    int p;
    pixel * myarr;
    pixel * myarr2;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    printf("Nprocs:\t%d, my_id:\t%d\n", p, my_id);

    if (my_id == ROOT)
    {
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

    /* Allocate a processing array for every node */
    myarr = (pixel *) malloc (elems_per_node * sizeof(pixel));

    /* Call MPI_Scatter to compute the sum of all pixels */
    MPI_Scatter((void *)src, elems_per_node, pixel_t_mpi, (void *)myarr, elems_per_node, pixel_t_mpi, ROOT, MPI_COMM_WORLD);
    uint mysum = get_px_sum(myarr, elems_per_node);

 /* We could have used MPI_AllReduce to avoid an extra communication step when
  * sending `avg` to the rest of the cores */
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
    {
        printf("Filtering took: %f secs\n", endtime - starttime);
        /* write result */
        printf("[ROOT] Writing output file\n");

        if (write_ppm (argv[2], xsize, ysize, (char *)src) != 0)
            exit(1);
    }

    /* Free allocated buffers */
    free(myarr);
    free(myarr2);
    free(src);

    /* Finalize MPI running environment */
    MPI_Finalize();

    return 0;
}
