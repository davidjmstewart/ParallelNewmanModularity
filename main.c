#include <stdio.h>
#include <stdlib.h>

#include "./lib/CDUTils.h"

unsigned long const MATRIX_SIZE = 500;
int main(int argc, char *argv[]) {
    int *A = (int *)malloc(MATRIX_SIZE * MATRIX_SIZE * sizeof(int));

    printf("Generating %lu x %lu adjacency matrix \n", MATRIX_SIZE, MATRIX_SIZE);
    genAdjacenyMatrix(A, MATRIX_SIZE);
}