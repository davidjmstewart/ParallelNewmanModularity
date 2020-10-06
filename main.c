#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "./lib/CDUTils.h"

unsigned long const MATRIX_SIZE = 20;
const int NUM_THREADS = 16;
long double elapsed;
struct timespec start, finish;

int main(int argc, char *argv[])
{
    int *A  = (int *)malloc(MATRIX_SIZE * MATRIX_SIZE * sizeof(int));   // our initial adjacency matrix that describes the graph
    double *B = (double *)malloc(MATRIX_SIZE * MATRIX_SIZE * sizeof(double)); // Modularity matrix
    // int A[25] = {0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0};
    // int A[25] = {0, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, -1, -1, 0};
    int *D = (int *)malloc(MATRIX_SIZE * sizeof(int));

    if (A == NULL) 
    {
        printf("Could not allocate %lu bytes to A. Exiting program\n", MATRIX_SIZE * MATRIX_SIZE * sizeof(int));
        return -1;
    }

    if (B == NULL)
    {
        printf("Could not allocate %lu bytes to B. Exiting program\n", MATRIX_SIZE * MATRIX_SIZE * sizeof(double));
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

    printf("Creating modularity Matrix \n");
    

    createModularityMatrix(B, A, D, MATRIX_SIZE, NUM_THREADS);
    printf("Outputting matlab file\n");
    matrixToFile(B, MATRIX_SIZE, MATLAB);
    // matrixToFile(B, MATRIX_SIZE, Python);



    printf("Performing power iteration \n");
    clock_gettime(CLOCK_MONOTONIC, &start);
    double *eigenVector = powerIteration(B, MATRIX_SIZE, NUM_THREADS, 0.001, 5000);
    // powerIteration(B, MATRIX_SIZE, NUM_THREADS, 50);
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("elapsed: %.9Lf \n", elapsed);
    printf("%0.9f\n", eigenVector[0]);
    membershipVectorToFile(eigenVector, MATRIX_SIZE, MATLAB);
    // membershipVectorToFile(eigenVector, MATRIX_SIZE, Python);
    // for (int i = 0; i < MATRIX_SIZE; i++)
    // {
    //     printf("%0.9f\n", eigenVector[i]);
    // }
    // free(A);
    free(D);
    free(B);

    return 0;
}

