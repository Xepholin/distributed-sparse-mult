#include "kernel.h"

#include <math.h>
#include <stdlib.h>

double norm_vec(const int n, const double a[n]) {
    double sum = 0;

    for (int i = 0; i < n; ++i) {
        sum += a[i] * a[i];
    }

    return sqrt(sum);
}

void dgemv(const CSR *a, const double *x, double *b) {
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

void pow_iter(const CSR *a, double *b_k, const int max_iter) {
    const int size = a->n;
    double *b_k1 = (double *)malloc(size * sizeof(double));

    for (int i = 0; i < max_iter; ++i) {
        dgemv(a, b_k, b_k1);
        const double b_k1_norm_inv = 1.0 / norm_vec(size, b_k1);

        for (int j = 0; j < size; ++j) {
            b_k[j] = b_k1[j] * b_k1_norm_inv;
        }
    }

    free(b_k1);
}