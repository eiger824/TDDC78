#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <mpi.h>

#include "ppmio.h"
#include "blurfilter-mpi.h"
#include "gaussw.h"
#include "image_utils.h"

#define MAX_RAD 1000
#define     ROOT    0

int main (int argc, char ** argv) {
    int radius;
    int xsize, ysize, colmax;
    pixel src[MAX_PIXELS];
    double w[MAX_RAD];

    /* Take care of the arguments */

    if (argc != 4)
    {
        fprintf(stderr, "Usage: %s radius infile outfile\n", argv[0]);
        exit(1);
    }

    radius = atoi(argv[1]);

    if (radius > MAX_RAD || radius < 1)
    {
        fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
        exit(1);
    }

    int my_id;
    int max_size;
    int p;
    int elems_per_node;
    pixel src_adj[MAX_PIXELS];
    pixel * myarr;
    pixel * myarr2;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    printf("Nprocs:\t%d, my_id:\t%d\n", p, my_id);

    if (my_id == ROOT)
    {
        /* read file */
        if (read_ppm (argv[2], &xsize, &ysize, &colmax, (char *) src) != 0)
            exit(1);

        if (colmax > 255)
        {
            fprintf(stderr, "Too large maximum color-component value\n");
            exit(1);
        }
        max_size = xsize * ysize;
        printf("[ROOT] Has read the image, generating coefficients\n");

        unsigned i;
        for (i = 1; i< p; ++i)
            MPI_Send(&max_size, 1, MPI_INT, i, ROOT, MPI_COMM_WORLD);
    }
    else
        MPI_Recv(&max_size, 1, MPI_INT, ROOT, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

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

    /* Gaussian weights */
    get_gauss_weights(radius, w);

    struct matrix original;
    struct matrix adjoint;

    double starttime, endtime;
    starttime = MPI_Wtime();

    /* Send xsize to the rest of the nodes */
    unsigned i;
    if (my_id == ROOT)
        for (i = 1; i < p; ++i)
            MPI_Send(&xsize, 1, MPI_INT, i, ROOT, MPI_COMM_WORLD);
    else
    {
        MPI_Recv(&xsize, 1, MPI_INT, ROOT, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        /* Update ysize now that we have both max_size and xsize */
        ysize = max_size / xsize;
    }

    //     elems_per_node = p * xsize;
    elems_per_node = (ysize / p) * xsize;

    if (my_id == ROOT)
    {
        printf("[ROOT] Every processing node with process %d elements.\n"
                , elems_per_node);
        printf("[ROOT] Calling filter on horizontal direction\n");
    }

    /* Synchronization point */
    MPI_Barrier(MPI_COMM_WORLD);

    /* Allocate a processing array for every node */
    myarr = (pixel * ) malloc (elems_per_node * sizeof(pixel));

    /* Scatter the work to do among processing units: horizontal filter first */
    MPI_Scatter((void *) src, elems_per_node, pixel_t_mpi,
            (void *) myarr, elems_per_node, pixel_t_mpi,
            ROOT, MPI_COMM_WORLD);

    blurfilter(myarr, elems_per_node, xsize, radius, w);
    /* Gather the results */
    MPI_Gather(myarr, elems_per_node, pixel_t_mpi,
            src, elems_per_node, pixel_t_mpi,
            ROOT, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    //     if (my_id == ROOT)
    //     {
    printf("[ROOT] Calculating adjoint pixel matrix\n");

    /* Convert our original image array `src` to a matrix */
    from_array(src, xsize, ysize, &original);

    /* Calculate its adjoint matrix */
    adj_matrix(&original, &adjoint);

    /* Convert the adjoint back to the array in `src_adj` */
    to_array(&adjoint, src_adj);

    printf("[ROOT] Calling filter on vertical direction.\n");
    //     }
    // The nodes must update their `elems_per_node` value
    elems_per_node = (xsize / p) * ysize;

    /* Allocate memory for the copy of the image */
    myarr2 = (pixel *) malloc (elems_per_node * sizeof(pixel));

    MPI_Scatter((void *) src_adj, elems_per_node, pixel_t_mpi,
            (void *) myarr2, elems_per_node, pixel_t_mpi,
            ROOT, MPI_COMM_WORLD);
    /* Call the filter */
    /* Note that xsize is now ysize since we have transposed the original matrix */
    blurfilter(src_adj, elems_per_node, ysize, radius, w);

    /* Gather the whole image into src_adj */
    MPI_Gather(myarr2, elems_per_node, pixel_t_mpi,
            src_adj, elems_per_node, pixel_t_mpi,
            ROOT, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    if (my_id == ROOT)
    {
        /* Finally, convert the adjoint array to the original array */
        from_array(src_adj, ysize, xsize, &adjoint);
        adj_matrix(&adjoint, &original);
        to_array(&original, src);
    }

    endtime = MPI_Wtime();

    if (my_id == ROOT)
    {
        printf("Filtering took: %f secs\n", endtime - starttime);

        /* write result */
        printf("Writing output file\n");

        if(write_ppm (argv[3], xsize, ysize, (char *)src) != 0)
            exit(1);

    }
    /* Free allocated buffers */
    free(myarr);
    free(myarr2);

    if (my_id == ROOT)
    {
        free_matrix(&original);
        free_matrix(&adjoint);
    }

    /* Finalize MPI running environment */
    MPI_Finalize();

    return 0;
}
