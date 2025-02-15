#define _POSIX_C_SOURCE 199309L

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
    int count = 0;
    int n = a->n;

    for (int i = 0; i < n; ++i) {
        double sum = 0.0;

        for (int j = a->row[count]; j < a->row[count + 1]; ++j) {
            sum += a->values[j] * x[i];
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

    // srand(time(NULL));
    srand(42);

    char *path = argv[1];
    int r = strtol(argv[2], NULL, 10);
    double elapsed = 0.0;
    struct timespec t1, t2;

    CSR *a = (CSR *)malloc(sizeof(CSR));
    double *x;
    double *b;

    import(a, path);

    x = (double *)malloc(a->n * sizeof(double));
    b = (double *)malloc(a->n * sizeof(double));

    init_matrix_r(a->n, x, 'r');
    init_matrix_r(a->n, b, 'z');

    clock_gettime(CLOCK_MONOTONIC_RAW, &t1);

    for (int i = 0; i < r; i++) {
        dgemv(a, x, b);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
    elapsed = (double)(t2.tv_nsec - t1.tv_nsec) / (double)r;

    for (int i = 0; i < a->n; ++i) {
        printf("%f ", b[i]);
    }

    printf("\n");
    printf("elpased: %lfs\n", elapsed * 1e-9);

    free(a->row);
    free(a->col);
    free(a->values);

    free(a);
    free(x);
    free(b);

    return 0;
}
