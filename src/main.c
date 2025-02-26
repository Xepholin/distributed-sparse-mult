#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "csr.h"
#include "distrib.h"
#include "kernel.h"
#include "utils.h"

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "./smult matrix.mtx iter\n");
        return 1;
    }

    MPI_Init(&argc, &argv);

    // srand(time(NULL));
    srand(42);

    char *path = argv[1];
    int r = strtol(argv[2], NULL, 10);
    int rank = 0;
    int n_ranks = 0;
    int size = 0;
    int nonzeros = 0;
    int last_rank = 0;

    CSR *a = (CSR *)malloc(sizeof(CSR));

    if (a == NULL) {
        perror("Error during allocations\n");
        exit(EXIT_FAILURE);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

    int n_rows[n_ranks];
    int offset_rows[n_ranks];
    int sendcounts[n_ranks];
    int offset[n_ranks];

    a->rank = rank;

    if (rank == 0) {
        printf("\nImporting matrix %s ...\n", path);
        import(a, path);
        size = a->n;
        nonzeros = a->total_values;

        last_rank = compute_partition(a, size, n_ranks, n_rows, offset_rows, sendcounts, offset);

        printf("Matrix %dx%d with %d values\n\n", size, size, a->total_values);
    }

    MPI_Bcast(&last_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Comm NEW_WORLD;

    if (rank <= last_rank) {
        MPI_Comm_split(MPI_COMM_WORLD, 1, rank, &NEW_WORLD);
    } else {
        MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, rank, &NEW_WORLD);
    }

    if (rank <= last_rank) {
        double *x;
        double *b_loc;
        double *b_glob;
        double *b;
        double *v;

        double lambda_loc = 0.0;
        double lambda_glob = 0.0;

        double mpi_elapsed_loc = 0.0;
        double mpi_elapsed = 0.0;
        double t1 = 0.0;
        double t2 = 0.0;

        int n = 0;

        MPI_Bcast(&size, 1, MPI_INT, 0, NEW_WORLD);
        MPI_Bcast(n_rows, n_ranks, MPI_INT, 0, NEW_WORLD);
        MPI_Bcast(offset_rows, n_ranks, MPI_INT, 0, NEW_WORLD);
        MPI_Bcast(sendcounts, n_ranks, MPI_INT, 0, NEW_WORLD);

        n = n_rows[rank] - 1;

        a->n = n;
        a->begin = offset_rows[rank];
        a->end = offset_rows[rank] + n;

        x = (double *)malloc(size * sizeof(double));
        b_glob = (double *)malloc(size * sizeof(double));
        b_loc = (double *)malloc(n * sizeof(double));
        b = (double *)malloc(size * sizeof(double));
        v = (double *)malloc(size * sizeof(double));

        if (x == NULL || b_glob == NULL || b_loc == NULL || b == NULL || v == NULL) {
            perror("Error during allocations\n");
            exit(EXIT_FAILURE);
        }

        init_matrix_r(size, x, 's');
        init_matrix_r(size, b_glob, 'z');
        init_matrix_r(n, b_loc, 'z');
        init_matrix_r(size, b, 'z');
        init_matrix_r(size, v, 'r');

        if (rank != 0) {
            a->row = (int *)malloc(n_rows[rank] * sizeof(int));

            if (a->row == NULL) {
                perror("Error during allocations\n");
                exit(EXIT_FAILURE);
            }
        }

        if (rank == 0) {
            printf("\nComputing ...\n");
        }

        MPI_Barrier(NEW_WORLD);
        t1 = MPI_Wtime();
        scatter_matrix_data(a, rank, n_rows, offset_rows, sendcounts, offset, NEW_WORLD);

        for (int i = 0; i < r; ++i) {
            dgemv(a, x, b_loc);
        }

        for (int i = 0; i < n_ranks; ++i) {
            n_rows[i]--;
        }

        MPI_Gatherv(b_loc, a->n, MPI_DOUBLE, b_glob, n_rows, offset_rows, MPI_DOUBLE, 0, NEW_WORLD);

        t2 = MPI_Wtime();
        mpi_elapsed_loc = t2 - t1;

        MPI_Reduce(&mpi_elapsed_loc, &mpi_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, NEW_WORLD);
        
        double v_norm = norm_vec(size, v);
        if (v_norm == 0.0) {
            perror("v norm too low\n");
        } else {
            double v_norm_inv = 1.0 / v_norm;
            for (int i = 0; i < size; ++i) {
                v[i] *= v_norm_inv;
            }

            lambda_loc = pow_iter(a, v, 1e-6, r);
            MPI_Reduce(&lambda_loc, &lambda_glob, 1, MPI_DOUBLE, MPI_MAX, 0, NEW_WORLD);
        }

        if (rank == 0) {
            double seq_elapsed = 0.0;

            a->n = size;
            a->total_values = nonzeros;

            t1 = MPI_Wtime();
            for (int i = 0; i < r; i++) {
                dgemv(a, x, b);
            }
            t2 = MPI_Wtime();
            seq_elapsed = t2 - t1;

            if (!check_results(size, b_glob, b)) {
                mpi_elapsed = (mpi_elapsed / (double)r);
                seq_elapsed = (seq_elapsed / (double)r);

                if (mpi_elapsed > 1e-4) {
                    printf("\nParallel elpased: %f s\n", mpi_elapsed);
                    printf("Sequential elpased: %f s\n", seq_elapsed);
                } else {
                    printf("\nParallel elpased: %f ns\n", mpi_elapsed * 1e9);
                    printf("Sequential elpased: %f ns\n", seq_elapsed * 1e9);
                }

                printf("Speedup: %f\n", (seq_elapsed / mpi_elapsed));
            }

            if (lambda_glob < 0) {
                printf("Error in Power Iteration\n");
            } else {
                printf("lambda value: %f\n", lambda_glob);
            }
        }

        free(a->row);
        free(a->col);
        free(a->values);

        free(x);
        free(b_loc);
        free(b_glob);
        free(b);
        free(v);

        MPI_Comm_free(&NEW_WORLD);
    }

    free(a);

    MPI_Finalize();

    return 0;
}
