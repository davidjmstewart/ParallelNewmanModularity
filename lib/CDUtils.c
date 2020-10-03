#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>

#include "./CDUTils.h"

// gcc -o CDUtils.o -c CDUtils.c

// Generates a size x size adjacency matrix in a flat structure
// INDEXING EXAMPLE
// If size is 5, the 3rd row, 4th column of A would be A[2*size + 3]
void genAdjacenyMatrix(int *A, unsigned long size) {
    srand(time(0));
    for (unsigned long i = 0; i < size; i++) {
        for (unsigned long j = i; j < size; j++) {
            if (i != j){ 
                int n = rand() % 2;
                A[i * size + j] = n;
                A[j * size + i] = n;
            } else
            { // diagonal must be 0
                A[i * size + j] = 0;
            }
        }
    }
}

void createDegreesVec(int *A, int *D, unsigned long size, int numThreads)
{
    unsigned long i, j;

    if (numThreads == USE_DEFAULT_NUM_THREADS) {
        numThreads = DEFAULT_NUM_THREADS;
    }

    omp_set_num_threads(numThreads);
    int sum = 0;
    #pragma omp parallel for private(i, j) reduction(+:sum)
    for (i = 0; i < size; i++)
    {
        int *ARow = &A[i * size]; // compute row once in this thread to save constant dereferencing
        for (j = 0; j < size; j++)
        {
            sum += *(ARow + j);
        }
        D[i] = sum;
    }
}

void createModularityMatrix(double *B, int *A, int *D, unsigned long size, int numThreads)
{
    printf("creating!\n");
    unsigned long i, j;
    int m = graphDegree(D, size)/2;

    if (numThreads == USE_DEFAULT_NUM_THREADS)
    {
        numThreads = DEFAULT_NUM_THREADS;
    }

    omp_set_num_threads(numThreads);
    #pragma omp parallel for private(i, j) 
    for (i = 0; i < size; i++)
    {
        // compute rows once in this thread to save constant dereferencing
        int *ARow = &A[i * size];
        double *BRow = &B[i * size];
        int iDegree = D[i];

        for (j = 0; j < size; j++)
        {
            BRow[j] = ARow[j] - (double)iDegree*D[j]/(2*m);
        }
    }
}

int graphDegree(int *D, unsigned long size) {
    int total = 0;
    for (int i = 0; i < size; i++)
    {
        total += D[i];
    }
    return total;
}
