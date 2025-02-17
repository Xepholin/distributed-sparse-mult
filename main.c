#define _POSIX_C_SOURCE 199309L

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef struct matrix_csr {
    int rank;
    int n;
    int total_values;
    int begin;
    int end;
    int *row;
    int *col;
    double *values;
} CSR;

void print_csr(const CSR *a) {
    printf("rank %d\n", a->rank);
    printf("row\n");

    for (int i = 0; i < a->n + 1; ++i) {
        printf("%d ", a->row[i]);
    }

    printf("\n\n");

    printf("col\n");
    for (int i = 0; i < a->total_values; ++i) {
        printf("%d ", a->col[i]);
    }

    printf("\n\n");

    printf("values\n");
    for (int i = 0; i < a->total_values; ++i) {
        printf("%.1lf ", a->values[i]);
    }

    printf("\n\n");
}

void init_matrix_r(int n, double a[n], char m) {
    // Random value per entry
    if (m == 'r' || m == 'R') {
        for (int i = 0; i < n; i++) {
            a[i] = (double)RAND_MAX / (double)rand();
        }
    } else {  // Zeroing up the array
        if (m == 'z' || m == 'Z') {
            for (int i = 0; i < n; i++) {
                a[i] = 0.0;
            }
        } else {  // Same value per entry
            if (m == 'c' || m == 'C') {
                double c = (double)RAND_MAX / (double)rand();
                for (int i = 0; i < n; i++) {
                    a[i] = c;
                }
            }
        }
    }
}

void import(CSR *a, const char *path) {
    FILE *file = fopen(path, "r");
    if (!file) {
        perror("Failed opening file");
        exit(EXIT_FAILURE);
    }

    char line[1024];

    while (fgets(line, sizeof(line), file)) {
        if (line[0] != '%') break;
    }

    int rows, cols, nonzeros;
    if (sscanf(line, "%d %d %d", &rows, &cols, &nonzeros) != 3) {
        perror("Wrong format need .mtx\n");
        exit(EXIT_FAILURE);
    }

    a->n = rows;
    a->total_values = nonzeros;
    a->values = (double *)malloc(nonzeros * sizeof(double));
    a->row = (int *)malloc((rows + 1) * sizeof(int));
    a->col = (int *)malloc(nonzeros * sizeof(int));

    int row = 0;
    int col = 0;
    int old_row = 0;
    int row_index = 0;
    double value = 0;

    a->row[row_index] = 0;

    for (int i = 0; i < nonzeros; ++i) {
        fscanf(file, "%d %d %lf\n", &col, &row, &value);

        row = row - 1;
        col = col - 1;

        if (row != old_row) {
            row_index++;
            a->row[row_index] = i;
            old_row = row;
        }

        a->col[i] = col;
        a->values[i] = value;
    }

    a->row[row_index + 1] = nonzeros;

    fclose(file);
}

void dgemv(const CSR *a, double *x, double *b) {
    int begin = a->begin;
    int count = 0;
    int n = a->n;

    for (int i = 0; i < n; ++i) {
        double sum = 0.0;

        for (int j = a->row[count]; j < a->row[count + 1]; ++j) {
            sum += a->values[j] * x[begin + i];
        }

        b[i] = sum;
        count++;
    }
}

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
        double elapsed_glob = 0.0;
        double seq_elapsed = 0.0;
        struct timespec t1, t2;

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

        init_matrix_r(size, x, 'r');
        init_matrix_r(size, b_glob, 'z');
        init_matrix_r(n, b_loc, 'z');
        init_matrix_r(size, b, 'z');

        if (rank != 0) {
            a->row = (int *)malloc(n_rows[rank] * sizeof(int));
        }

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

        for (int i = 0; i < r; ++i) {
            clock_gettime(CLOCK_MONOTONIC_RAW, &t1);

            dgemv(a, x, b_loc);

            clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
            MPI_Barrier(NEW_WORLD);

            mpi_elapsed_loc += ((double)t2.tv_sec - (double)t1.tv_sec) + ((double)t2.tv_nsec - (double)t1.tv_nsec);
        }

        for (int i = 0; i < n_ranks; ++i) {
            n_rows[i]--;
        }

        MPI_Gatherv(b_loc, a->n, MPI_DOUBLE, b_glob, n_rows, offset_rows, MPI_DOUBLE, 0, NEW_WORLD);
        MPI_Reduce(&mpi_elapsed_loc, &elapsed_glob, 1, MPI_DOUBLE, MPI_SUM, 0, NEW_WORLD);

        if (rank == 0) {
            a->n = size;
            a->total_values = nonzeros;

            for (int i = 0; i < r; i++) {
                clock_gettime(CLOCK_MONOTONIC_RAW, &t1);

                dgemv(a, x, b);

                clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
                seq_elapsed += ((double)t2.tv_sec - (double)t1.tv_sec) + ((double)t2.tv_nsec - (double)t1.tv_nsec);
            }
        }

        if (rank == 0) {
            int count = 0;
            for (int i = 0; i < size; ++i) {
                if (b_glob[i] == b[i]) {
                    count++;
                } else {
                    printf("%f %f\n", b_glob[i], b[i]);
                }
            }

            printf("\nValidation DGEMV...\n");

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

            printf("\nParallel elpased: %lf s\n", (elapsed_glob / (double)r) * 10e-9);
            printf("Sequential elpased: %lf s\n", (seq_elapsed / (double)r) * 10e-9);
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
