#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "./CDUTils.h"

// gcc -o CDUtils.o -c CDUtils.c

// Generates a size x size adjacency matrix in a flat structure
// INDEXING EXAMPLE
// If size is 5, the 3rd row, 4th column of A would be A[2*size + 3]
void genAdjacenyMatrix(int *A, unsigned long size) {
    for (unsigned long i = 0; i < size; i++) {
        for (unsigned long j = i; j < size; j++) {
            int n = rand() % 2;
            A[i * size + j] = n;
            A[j * size + i] = n;
        }
    }
}

void createDegreesVec(int *A, int *D, unsigned long size, int numThreads)
{
    unsigned long i, j, k;
    double elapsed;

    omp_set_num_threads(numThreads);
    int sum;
#pragma omp parallel for private(i, k, tmp) reduction(+:sum)
    for (i = 0; i < size; i++)
    {
        sum = 0;
        int *ARow = &A[i * size]; // compute row once in this thread to save constant dereferencing
        for (k = 0; k < size; k++)
        {
            sum += *(ARow + k);
        }
        D[i] = sum;
    }
}
