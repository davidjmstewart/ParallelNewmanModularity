#include <time.h>
#include <stdio.h>
#include <math.h>

#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <omp.h>
#include <assert.h>

#define THREAD_RANGE 16 // Run for 1:THREAD_RANGE threads
#define NUM_AVERAGES 10 // take the average of 5 timings for each matrix size, and each number of threads
#define MATRIX_SIZE  3000 

// gcc -fopenmp matrix_tasking.c -o matrix_tasking -O3 -Wall -Werror
double parallelTimings[THREAD_RANGE];

void matrixVectorMultiplicationSequential(double *restrict M, double *restrict V, double *restrict results,  unsigned long matrixSize)
{
    int i, j;
    for (i = 0; i < matrixSize; i++)
    {
        double *MHead = &M[i * matrixSize];
        double tmp = 0;
        for (j = 0; j < matrixSize; j++)
        {
            tmp += MHead[j] * V[j];
        }
        results[i] = tmp;
    }
}


void matrixVectorMultiplicationParallel(double *restrict M, double *restrict V, double *restrict results, unsigned long matrixSize, int numThreads)
{
    omp_set_num_threads(numThreads);
    unsigned long i, j;

    #pragma omp parallel for private(j)
    for (i = 0; i < matrixSize; i++)
    {
        double tmp = 0;
        double *MHead = &M[i * matrixSize];

        #pragma omp simd reduction(+ \
                          : tmp)
        for (j = 0; j < matrixSize; j++)
        {
            tmp += MHead[j] * V[j];
        }
        results[i] = tmp; 
    }
}

void doParallelComputation(double *restrict A, double *restrict B, double *restrict V, double *restrict resultsA, double *restrict resultsB, unsigned long matrixSize, int numThreads)
{
    #pragma omp task
    matrixVectorMultiplicationParallel(A,  V, resultsA, matrixSize, numThreads/2);
    #pragma omp task
    matrixVectorMultiplicationParallel(B, V, resultsB, matrixSize, numThreads/2);
    #pragma omp taskwait
}

void genRandVector(double *S, unsigned long size)
{
    srand(time(0));
    unsigned long i;
    #pragma omp parallel for private(i)
    for (i = 0; i < size; i++)
    {
        double n = rand() % 3;
        S[i] = n;
    }
}

void doSequentialComputation(double *restrict A, double *restrict B, double *restrict V, double *restrict resultsA, double *restrict resultsB, unsigned long matrixSize)
{
    matrixVectorMultiplicationSequential(A, V, resultsA, matrixSize);
    matrixVectorMultiplicationSequential(B, V, resultsB, matrixSize);
}


void genRandMatrix(double *A, unsigned long size)
{
    srand(time(0));
    unsigned long i, j;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            double n = rand() % 3;
            A[i * size + j] = n;
        }
    }
}

int main(int argc, char *argv[])
{
    struct timespec start, finish;
    double elapsed;

    unsigned long matrixSize = 100;

    double *V = (double *)malloc(matrixSize * sizeof(double));
    double *seqVA = (double *)malloc(matrixSize * sizeof(double)); // Store the results of A*v in the sequential implementation here
    double *parVA = (double *)malloc(matrixSize * sizeof(double)); // Store the results of A*v in the parallel implementation here
    double *seqVB = (double *)malloc(matrixSize * sizeof(double)); // Store the results of B*v in the sequential implementation here
    double *parVB = (double *)malloc(matrixSize * sizeof(double)); // Store the results of B*v in the parallel implementation here
    double *A = (double *)malloc(matrixSize * matrixSize * sizeof(double)); // First matrix to multiply by V
    double *B = (double *)malloc(matrixSize * matrixSize * sizeof(double)); // Second matrix to multiply by V

    genRandVector(V, matrixSize);
    genRandMatrix(A, matrixSize);
    genRandMatrix(B, matrixSize);
    double sequentialTiming = 0;
    
    for (int a = 0; a < NUM_AVERAGES; a++)
    {
        clock_gettime(CLOCK_MONOTONIC, &start);
        doSequentialComputation(A, B, V, seqVA, seqVB, matrixSize);
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        sequentialTiming += elapsed;
        for (int t = 1; t <= THREAD_RANGE; t++)
        {
            clock_gettime(CLOCK_MONOTONIC, &start);
            omp_set_num_threads(t);

            #pragma omp parallel 
            {
                #pragma omp single 
                doParallelComputation(A, B, V, parVA, parVB, matrixSize, t);
            }
                // doParallelComputation(A, B, V, parVA, parVB, matrixSize, t);

            clock_gettime(CLOCK_MONOTONIC, &finish);
            elapsed = (finish.tv_sec - start.tv_sec);
            elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
            parallelTimings[t - 1] += elapsed;

            // parallelTiming += elapsed;
            for (int i = 0; i < matrixSize; i++)
            {
                assert(fabs(seqVA[i] - parVA[i]) < 0.01);
            }
            
        }
    }
    sequentialTiming /= NUM_AVERAGES;

    printf("Sequential: %f \n", sequentialTiming);
    printf("Parallel: ");

    for (int t = 0; t < THREAD_RANGE; t++)
    {
        parallelTimings[t] /= NUM_AVERAGES;
        printf("%f ", parallelTimings[t]);
    }
    printf("\n");

    free(seqVA);
    free(parVA);
    free(seqVB);
    free(parVB);
    free(A);
    free(V);
    free(B);

    return 0;
}