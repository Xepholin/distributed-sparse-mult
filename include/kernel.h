#pragma once

#include "csr.h"

double norm_vec(const int n, const double a[n]);
void dgemv(const CSR *a, const double *x, double *b);
double pow_iter(const CSR *a, double *v_k, const double threshold, const int max_iter);