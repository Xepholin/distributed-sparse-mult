#pragma once

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

void import(CSR *a, const char *path);