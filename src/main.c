#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "csr.h"
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

        int n_mod = size % n_ranks;
        int n_div = size / n_ranks;

        int cumul = 0;

        for (int p = 0; p < n_ranks; ++p) {
            n_rows[p] = n_div + (p < n_mod ? 1 : 0);
            offset_rows[p] = cumul;
            cumul += n_rows[p];

            n_rows[p]++;

            if (offset_rows[p] + n_rows[p] >= size + 1) {
                last_rank = p;
                break;
            }
        }

        int index = 0;
        for (int p = 0; p < n_ranks; ++p) {
            int count = a->row[offset_rows[p] + n_rows[p] - 1] - a->row[offset_rows[p]];

            sendcounts[p] = count;
            offset[p] = index;

            index += count;

            if (sendcounts[p] + offset[p] >= a->row[size]) {
                break;
            }
        }

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

        int n = 0;

        double mpi_elapsed_loc = 0.0;
        double mpi_elapsed = 0.0;
        double seq_elapsed = 0.0;
        double t1 = 0.0;
        double t2 = 0.0;

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

        init_matrix_r(size, x, 's');
        init_matrix_r(size, b_glob, 'z');
        init_matrix_r(n, b_loc, 'z');
        init_matrix_r(size, b, 'z');

        if (rank != 0) {
            a->row = (int *)malloc(n_rows[rank] * sizeof(int));
        }

        t1 = MPI_Wtime();
        MPI_Scatterv(a->row, n_rows, offset_rows, MPI_INT, a->row, n_rows[rank], MPI_INT, 0, NEW_WORLD);

        a->total_values = a->row[a->n] - a->row[0];

        printf("Process %d: processing of %d values over %d rows.\n", rank, a->total_values, a->n);

        if (rank != 0) {
            int temp = a->row[0];

            for (int i = 0; i < a->n + 1; ++i) {
                a->row[i] -= temp;
            }

            a->col = (int *)malloc(a->total_values * sizeof(int));
            a->values = (double *)malloc(a->total_values * sizeof(double));
        }

        MPI_Scatterv(a->col, sendcounts, offset, MPI_INT, a->col, sendcounts[rank], MPI_INT, 0, NEW_WORLD);
        MPI_Scatterv(a->values, sendcounts, offset, MPI_DOUBLE, a->values, sendcounts[rank], MPI_DOUBLE, 0, NEW_WORLD);

        if (rank == 0) {
            printf("\nComputing ...\n");
        }

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
        MPI_Barrier(NEW_WORLD);

        if (rank == 0) {
            a->n = size;
            a->total_values = nonzeros;

            t1 = MPI_Wtime();
            for (int i = 0; i < r; i++) {
                dgemv(a, x, b);
            }
            t2 = MPI_Wtime();
            seq_elapsed = t2 - t1;

            printf("\nChecking result...\n");
            int count = 0;
            for (int i = 0; i < size; ++i) {
                if (b_glob[i] == b[i]) {
                    count++;
                } else {
                    printf("%f %f\n", b_glob[i], b[i]);
                }
            }

            if (count == size) {
                printf("DGEMV passed.\n");
            } else {
                printf("DGEMV failed.\n");
                printf("Number of wrong element: %d/%d\n", count, size);

                for (int i = 0; i < size; ++i) {
                    if (b_glob[i] == b[i]) {
                        continue;
                    } else {
                        printf("Element %d: %f != %f\n", i, b_glob[i], b[i]);
                    }
                }
            }

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

        free(a->row);
        free(a->col);
        free(a->values);

        free(x);
        free(b_loc);
        free(b_glob);
        free(b);

        MPI_Comm_free(&NEW_WORLD);
    }

    free(a);

    MPI_Finalize();

    return 0;
}
