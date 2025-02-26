#pragma once

#include "csr.h"

double norm_vec(const int n, const double a[n]);
void dgemv(const CSR *a, const double *x, double *b);
void pow_iter(const CSR *a, double *b_k, const int max_iter);