#include <stdio.h>
#include <stdlib.h>

#include "./lib/CDUTils.h"

unsigned long const MATRIX_SIZE = 5000;

int main(int argc, char *argv[]) {
    int *A = (int *)malloc(MATRIX_SIZE * MATRIX_SIZE * sizeof(int));
    if (A == NULL) {
        printf("Could not allocate %lu bytes to A. Exiting program\n", MATRIX_SIZE * MATRIX_SIZE * sizeof(int));
        return -1;
    }
    printf("Generating %lu x %lu adjacency matrix \n", MATRIX_SIZE, MATRIX_SIZE);
    genAdjacenyMatrix(A, MATRIX_SIZE);

    free(A);
    return 0;
}