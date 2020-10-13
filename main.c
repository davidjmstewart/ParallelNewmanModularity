#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "./lib/CDUTils.h"

unsigned long const MATRIX_SIZE = 10;
const int NUM_THREADS = 2;
long double elapsed;
struct timespec start, finish;

int main(int argc, char *argv[])
{
    int *A  = (int *)malloc(MATRIX_SIZE * MATRIX_SIZE * sizeof(int));   // our initial adjacency matrix that describes the graph
    double *B = (double *)malloc(MATRIX_SIZE * MATRIX_SIZE * sizeof(double)); // Modularity matrix

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
    genAdjacencyMatrix(A, MATRIX_SIZE);

    printf("Creating degree vector \n");
    createDegreesVec(A, D, MATRIX_SIZE, NUM_THREADS);

    printf("Creating modularity Matrix \n");
    createModularityMatrix(B, A, D, MATRIX_SIZE, NUM_THREADS);
    printf("Outputting matlab file\n");
    matrixToFile(B, MATRIX_SIZE, MATLAB);
    // matrixToFile(B, MATRIX_SIZE, Python);
    integerMatrixToFile(A, MATRIX_SIZE, MATLAB);
    printf("Performing power iteration \n");
    clock_gettime(CLOCK_MONOTONIC, &start);
    // eigenPair eigenPair = powerIteration(B, MATRIX_SIZE, NUM_THREADS, 0.0000001, 5000);
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("elapsed: %.9Lf \n", elapsed);
    // printf("%0.9f\n", eigenPair.eigenvalue);
    // membershipVectorToFile(eigenPair.eigenvector, MATRIX_SIZE, MATLAB);
    // membershipVectorToFile(eigenVector, MATRIX_SIZE, Python);
    int *globalVertices = (int *)malloc(MATRIX_SIZE * sizeof(int)); // Store the subgraph modularity matrix computed in parallel
    createGlobalVertices(globalVertices, MATRIX_SIZE, 4);
    assignCommunity(B, MATRIX_SIZE, MATRIX_SIZE, globalVertices, NUM_THREADS);
    free(A);
    free(D);
    free(B);
    free(globalVertices);

    return 0;
}

