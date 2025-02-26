#include "csr.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

    if (a->values == NULL || a->row == NULL || a->col == NULL) {
        perror("Error during allocations\n");
        exit(EXIT_FAILURE);
    }

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