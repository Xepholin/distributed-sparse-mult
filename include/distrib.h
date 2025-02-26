#pragma once

#include <mpi.h>

#include "csr.h"

int compute_partition(const CSR *a, int size, int n_ranks, int n_rows[n_ranks],
                      int offset_rows[n_ranks], int sendcounts[n_ranks],
                      int offset[n_ranks]);

void scatter_matrix_data(CSR *a, int rank, int *n_rows, int *offset_rows, int *sendcounts, int *offset, MPI_Comm NEW_WORLD);

int check_results(int size, double *b_glob, double *b);