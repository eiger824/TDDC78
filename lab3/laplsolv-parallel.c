#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

double find_max(double * array, size_t n)
{
    double max = array[0];
    for (size_t i = 1; i < n; ++i)
    {
        if (array[i] > max)
            max = array[i];
    }
    return max;
}

void init_big_T(double ** T, int rows, int cols)
{
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            // Always, j=0 or j==n-1 means a 1 must be set
            if (j == 0 || j == cols - 1)
            {
                if (i < rows - 1) 
                    T[i][j] = 1.0;
                else
                    T[i][j] = 2.0;
            }
            else
            {
                if (i == rows - 1)
                    T[i][j] = 2.0;
                else
                    T[i][j] = 0.0;
            }
        }
    }
}

int main(int argc, char* argv[])
{
    int nr_threads;
    int n = 1000;
    int j;
    int my_id;
    int maxiter = 1000;
    double tol = 1.0e-3;
    /* Allocate a big matrix */
    double ** T = (double **) malloc (sizeof *T * (n + 2));
    for (int i = 0 ; i < n + 2; ++i)
        *(T + i) = (double *) malloc (sizeof(double) * (n + 2));
    double error, t0, t1;
    int * start_pos;
    int * end_pos;
    double *tmp2;
#pragma omp parallel
    {
        nr_threads = omp_get_num_threads();

        /* Allocate tmp2 - array */
        tmp2 = (double *) malloc (sizeof *tmp2 * n);
    }
    /* Allocate tmp's 1 and 3 - matrices */
    double ** tmp1 = (double ** ) malloc (sizeof *tmp1 * n);
    double ** tmp3 = (double ** ) malloc (sizeof *tmp3 * n);
    for (int i = 0; i < n; ++i)
        *(tmp1 + i) = (double *) malloc (sizeof (double) * nr_threads);
    for (int i = 0; i < n; ++i)
        *(tmp3 + i) = (double *) malloc (sizeof (double) * nr_threads);

    /* Set boundary conditions and initial values for the unknowns */
    init_big_T(T, n + 2, n + 2);

    int ratio = n / nr_threads;

    // Init the start & end positions
    start_pos = (int * ) malloc (sizeof *start_pos * nr_threads);
    end_pos   = (int * ) malloc (sizeof *end_pos * nr_threads);

    for (int i = 0; i < nr_threads; ++i)
    {
        start_pos[i] = (i * ratio) + 1;
        end_pos[i]   = (i + 1) * ratio;
    }

    t0 = omp_get_wtime();

    int iter; 
    for (iter = 0; iter <= maxiter; ++iter)
    {
        error = 0.0;

        for (int i = 0; i < nr_threads; ++i) // cols in big T
        {
            for (int k = 1; k <= n; ++k) // rows in big T
            {
                tmp1[k-1][i] = T[k][start_pos[i] - 1];
                tmp3[k-1][i] = T[k][end_pos[i] + 1];
            }
        }
#pragma omp parallel private(j,my_id) shared(T,tmp1,tmp3,tmp2)
        {
            my_id = omp_get_thread_num();
            double * diff = (double *) malloc(sizeof *diff * n);
            double current_error = 0.0;
            for ( int pos = start_pos[my_id]; pos < end_pos[my_id]; ++pos )
            {
                for (int k = 0; k < n; ++k)
                {
                    tmp2[k] = T[k+1][pos];
                    T[k+1][pos] = (T[k][pos] + T[k+2][pos] + T[k+1][pos+1] + tmp1[k][my_id]) / 4.0;
                    diff[k]  = fabs(tmp2[k] - T[k+1][pos]);
                    tmp1[k][my_id] = tmp2[k];
                }
                current_error = find_max(diff, n);
                /* Get the maximum */
                error = error > current_error ? error : current_error;
            }
            for ( int k = 0; k < n; ++k)
            {
                tmp2[k] = T[k+1][end_pos[my_id]];
                T[k+1][end_pos[my_id]] =
                    (T[k][end_pos[my_id]] + T[k+2][end_pos[my_id]] + tmp3[k][my_id] + tmp1[k][my_id]) / 4.0;
                diff[k] = fabs(tmp2[k] - T[k+1][end_pos[my_id]]);
                tmp1[k][my_id] = tmp2[k];
            }
            current_error = find_max(diff, n);
            error = error > current_error ? error : current_error;

            free(diff);
        }
        if (error < tol)
            break;
    }

    t1 = omp_get_wtime();

    printf("-----------------------\n");
    printf("Elapsed:\t%.6f s.\n", t1 - t0);
    printf("Iterations:\t%d\n", iter);
    printf("Temperature:\t%.12f\n", T[1][1]);
    printf("Nr. threads:\t%d\n", nr_threads);
    printf("-----------------------\n");

    /* Free stuff */
    for (int i = 0; i < n + 2; ++i)
        free(*(T + i));
    free(T);

    return 0;
}

