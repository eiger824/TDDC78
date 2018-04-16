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

int main (int argc, char ** argv) {
    int radius;
    int xsize, ysize, colmax;
    pixel src[MAX_PIXELS];
    struct timespec stime, etime;
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

#ifdef WITH_MPI
    int p;
    int my_id;
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

    printf("Calling filter\n");

    /* Different ways of measuring time, according to which technique to use */
#ifdef WITH_MPI
    double starttime, endtime;
    starttime = MPI_Wtime();
#endif

#ifdef WITH_OPENMP
    clock_gettime(CLOCK_REALTIME, &stime);
#endif

    /* Apply the filter */
    blurfilter(xsize, ysize, src, radius, w);

#ifdef WITH_MPI
    endtime = MPI_Wtime();
    printf("Filtering took: %f secs\n", endtime - starttime);
#endif

#ifdef WITH_OPENMP
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
