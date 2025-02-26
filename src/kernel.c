#include "kernel.h"
#include "utils.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double norm_vec(const int n, const double a[n]) {
    double max = 0.0;

    for (int i = 0; i < n; ++i) {
        if (fabs(a[i]) > max) {
            max = fabs(a[i]);
        }
    }
    return max;
}

void dgemv(const CSR *a, const double *x, double *b) {
    const int begin = a->begin;
    const int n = a->n;
    int count = 0;

    for (int i = 0; i < n; ++i) {
        double sum = 0.0;

        for (int j = a->row[count]; j < a->row[count + 1]; ++j) {
            sum += a->values[j] * x[begin + i];
        }

        b[i] = sum;
        count++;
    }
}

double pow_iter(const CSR *a, double *v_k, const double threshold, const int max_iter) {
    const int n = a->n;
    const int begin = a->begin;
    double *v_k1 = (double *)malloc(n * sizeof(double));

    if (v_k1 == NULL) {
        return -1;
    }

    init_matrix_r(n, v_k1, 'z');

    double lambda_prev = 0.0;

    for (int i = 0; i < max_iter; ++i) {
        dgemv(a, v_k, v_k1);

        double lambda = norm_vec(n, v_k1);

        double v_k1_norm_inv = 1.0 / lambda;
        for (int j = 0; j < n; ++j) {
            v_k[begin + j] = v_k1[j] * v_k1_norm_inv;
        }

        if (fabs(lambda - lambda_prev) < threshold) {
            free(v_k1);
            return lambda;
        }

        lambda_prev = lambda;
    }

    free(v_k1);
    return lambda_prev;
}