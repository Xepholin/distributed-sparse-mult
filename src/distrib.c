#include "distrib.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "kernel.h"
#include "utils.h"

#define TOLERANCE 1e-12

int compute_partition(const CSR *a, int size, int n_ranks, int n_rows[n_ranks],
                      int offset_rows[n_ranks], int sendcounts[n_ranks],
                      int offset[n_ranks]) {
    int n_mod = size % n_ranks;
    int n_div = size / n_ranks;
    int cumul = 0;
    int last_rank = 0;

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

    return last_rank;
}

void scatter_matrix_data(CSR *a, int rank, int *n_rows, int *offset_rows, int *sendcounts, int *offset, MPI_Comm NEW_WORLD) {
    MPI_Scatterv(a->row, n_rows, offset_rows, MPI_INT, a->row, n_rows[rank], MPI_INT, 0, NEW_WORLD);
    a->total_values = a->row[a->n] - a->row[0];

    if (rank != 0) {
        int temp = a->row[0];
        for (int i = 0; i < a->n + 1; ++i) {
            a->row[i] -= temp;
        }
        a->col = (int *)malloc(a->total_values * sizeof(int));
        a->values = (double *)malloc(a->total_values * sizeof(double));
    }

    printf("Process %d: processing of %d values over %d rows.\n", rank, a->total_values, a->n);

    MPI_Scatterv(a->col, sendcounts, offset, MPI_INT, a->col, sendcounts[rank], MPI_INT, 0, NEW_WORLD);
    MPI_Scatterv(a->values, sendcounts, offset, MPI_DOUBLE, a->values, sendcounts[rank], MPI_DOUBLE, 0, NEW_WORLD);
}

int check_results(int size, double *b_glob, double *b) {
    printf("\nChecking result...\n");
    int count = 0;

    for (int i = 0; i < size; ++i) {
        if (fabs(b_glob[i] - b[i]) < TOLERANCE) {
            count++;
        } else {
            printf("Mismatch at element %d: %f != %f\n", i, b_glob[i], b[i]);
        }
    }

    if (count == size) {
        printf("DGEMV passed.\n");
        return 0;
    } else {
        printf("DGEMV failed.\n");
        printf("Number of correct elements: %d/%d\n", count, size);
        return 1;
    }
}
