#include "utils.h"

#include <stdio.h>
#include <stdlib.h>

void print_csr(const CSR *a) {
    printf("rank %d\n", a->rank);

    printf("n: %d\n", a->n);
    printf("nonzeros: %d\n", a->total_values);

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
                // double c = (double)RAND_MAX / (double)rand();
                double c = 1.0;
                for (int i = 0; i < n; i++) {
                    a[i] = c;
                }
            } else {
                if (m == 's' || m == 'S') {
                    for (int i = 0; i < n; ++i) {
                        a[i] = (double)i + 1.0;
                    }
                }
            }
        }
    }
}