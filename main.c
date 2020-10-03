#include <stdio.h>
#include <stdlib.h>

#include "./lib/CDUTils.h"

unsigned long const MATRIX_SIZE = 3;
const int NUM_THREADS = 16;

int main(int argc, char *argv[])
{
    int *A = (int *)malloc(MATRIX_SIZE * MATRIX_SIZE * sizeof(int));
    int *D = (int *)malloc(MATRIX_SIZE * sizeof(int));

    if (A == NULL) 
    {
        printf("Could not allocate %lu bytes to A. Exiting program\n", MATRIX_SIZE * MATRIX_SIZE * sizeof(int));
        return -1;
    }

    if (D == NULL)
    {
        printf("Could not allocate %lu bytes to D. Exiting program\n", MATRIX_SIZE * sizeof(int));
        return -1;
    }


    printf("Generating %lu x %lu adjacency matrix \n", MATRIX_SIZE, MATRIX_SIZE);
    genAdjacenyMatrix(A, MATRIX_SIZE);

    printf("Creating degree vector \n");

    createDegreesVec(A, D, MATRIX_SIZE, NUM_THREADS);

    free(A);
    free(D);

    return 0;
}