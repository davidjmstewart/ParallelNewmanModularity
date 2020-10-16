#include <time.h>
#include <stdio.h>
#include <math.h>

#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <omp.h>
#include <assert.h>

#define N 40000 // Matrix size will be N x N
#define T 1
#define THREAD_RANGE 16 // Run for 1:THREAD_RANGE threads
#define NUM_AVERAGES 10 // take the average of 5 timings for each matrix size, and each number of threads
#define NUM_MATRIX_SIZES 8
// unsigned long matrixSizes[NUM_MATRIX_SIZES] = {51200};

unsigned long matrixSizes[NUM_MATRIX_SIZES] = {100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200};
// unsigned long matrixSizes[NUM_MATRIX_SIZES] = {100000, 200000, 500000, 1000000, 2000000, 3000000, 5000000, 10000000, 20000000, 50000000};

// unsigned long matrixSizes[NUM_MATRIX_SIZES] = { 51200 };

double sequentialTimings[NUM_MATRIX_SIZES];
double parallelTimings[NUM_MATRIX_SIZES][THREAD_RANGE];

//gcc -fopenmp -g ./benchmarking/matrix_vector_multiplication/mvm_benchmarking.c -o ./benchmarking/matrix_vector_multiplication/mvm ./lib/CDUtils.o -lm -Wall -Wpedantic -Waggressive-loop-optimizations -O3 -march=native -fopt-info-vec-missed=v.txt
// #pragma GCC option("arch=native", "tune=native", "no-zero-upper") //Enable AVX
#pragma GCC target("avx") //Enable AVX


void matrixVectorMultiplicationSequential(double *restrict M, double *restrict V, double *restrict results,  unsigned long matrixSize)
{
    // double *V3 = (double *)malloc(2 * matrixSize * sizeof(double));
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

#pragma omp reduction(+ \
                           : tmp)
        for (j = 0; j < matrixSize; j++)
        {
            //results[i] += A[i * matrixSize + j] * V[j];
            tmp += MHead[j] * V[j];
        }
        results[i] = tmp; // write-only to results, not adding to old value.
    }
}

void doParallelComputation(double *restrict A, double *restrict B, double *restrict V, double *restrict resultsA, double *restrict resultsB, unsigned long matrixSize, int numThreads)
{

// #pragma omp task
    matrixVectorMultiplicationParallel(A,  V, resultsA, matrixSize, numThreads);
// #pragma omp task
    matrixVectorMultiplicationParallel(B, V, resultsB, matrixSize, numThreads);

// #pragma omp taskwait
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

    for (int m = 0; m < NUM_MATRIX_SIZES; m++)
    {
        sequentialTimings[m] = 0;
        for (int t = 0; t < THREAD_RANGE; t++)
        {
            parallelTimings[m][t] = 0;
        }
    }

    for (int m = 0; m < NUM_MATRIX_SIZES; m++)
    {

        unsigned long matrixSize = matrixSizes[m]; //

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

        for (int a = 0; a < NUM_AVERAGES; a++)
        {
            clock_gettime(CLOCK_MONOTONIC, &start);
            doSequentialComputation(A, B, V, seqVA, seqVB, matrixSize);
            clock_gettime(CLOCK_MONOTONIC, &finish);
            elapsed = (finish.tv_sec - start.tv_sec);
            elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
            sequentialTimings[m] += elapsed;
            for (int t = 1; t <= THREAD_RANGE; t++)
            {
                clock_gettime(CLOCK_MONOTONIC, &start);
                #pragma omp parallel 
                {
                #pragma omp single // we want a single thread to enter this initially
                    doParallelComputation(A, B, V, parVA, parVB, matrixSize, t);
                }

                clock_gettime(CLOCK_MONOTONIC, &finish);
                elapsed = (finish.tv_sec - start.tv_sec);
                elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
                parallelTimings[m][t - 1] += elapsed;
                for (int i = 0; i < matrixSize; i++)
                {
 
                        // printf("i: %d \t j: %d \n", i, j);
                        // printf("seqVA: %.19f \t parVA: %.19f \n", seqVA[i * matrixSize + j], parVA[i * matrixSize + j]);

                        assert(fabs(seqVA[i ] - parVA[i ]) < 0.01);
                        assert(fabs(seqVB[i ] - parVB[i ]) < 0.01);
                }
                
            }
        }
        free(seqVA);
        free(parVA);
        free(seqVB);
        free(parVB);
        free(A);
        free(V);
        free(B);

    }

    for (int m = 0; m < NUM_MATRIX_SIZES; m++)
    {
        sequentialTimings[m] /= NUM_AVERAGES;
        for (int t = 0; t < THREAD_RANGE; t++)
        {
            parallelTimings[m][t] /= NUM_AVERAGES;
        }
    }

    printf("Sequential:\n");
    for (int m = 0; m < NUM_MATRIX_SIZES; m++)
    {
        printf("%lux%lu:\t", matrixSizes[m], matrixSizes[m]);

        printf("%.9f \n", sequentialTimings[m]);
    }

    printf("\nParallel:\n");
    for (int m = 0; m < NUM_MATRIX_SIZES; m++)
    {
        printf("%lux%lu:\t", matrixSizes[m], matrixSizes[m]);

        for (int t = 0; t < THREAD_RANGE; t++)
        {
            printf(" %.9f ", parallelTimings[m][t]);
        }
        printf("\n");
    }
    printf("\n");

    return 0;
}