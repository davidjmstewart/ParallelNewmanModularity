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
#define NUM_MATRIX_SIZES 6
// unsigned long matrixSizes[NUM_MATRIX_SIZES] = {51200};

unsigned long matrixSizes[NUM_MATRIX_SIZES] = { 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200 };
// unsigned long matrixSizes[NUM_MATRIX_SIZES] = {100000, 200000, 500000, 1000000, 2000000, 3000000, 5000000, 10000000, 20000000, 50000000};

// unsigned long matrixSizes[NUM_MATRIX_SIZES] = { 51200 };

double sequentialTimings[NUM_MATRIX_SIZES];
double parallelTimings[NUM_MATRIX_SIZES][THREAD_RANGE];

//gcc -fopenmp -g ./benchmarking/matrix_vector_multiplication/mvm_benchmarking.c -o ./benchmarking/matrix_vector_multiplication/mvm ./lib/CDUtils.o -lm -Wall -Wpedantic -Waggressive-loop-optimizations -O3 -march=native -fopt-info-vec-missed=v.txt
// #pragma GCC option("arch=native", "tune=native", "no-zero-upper") //Enable AVX
#pragma GCC target("avx")                                         //Enable AVX

void doSequentialComputation(double *A, double *V, double *results, unsigned long matrixSize)
{
    // double *V3 = (double *)malloc(2 * matrixSize * sizeof(double));
    int i, j;

    for (i = 0; i < matrixSize; i++)
    {
        double *AHead = &A[i * matrixSize];
        double tmp = 0;

        for (j = 0; j < matrixSize; j++)
        {
            tmp += AHead[j] * V[j];
        }
        results[i] = tmp;
    }
}

void doParallelComputation(double *restrict A, double *restrict V, double *restrict results, unsigned long matrixSize, int numThreads)
{
    omp_set_num_threads
    const int BLOCK_SIZE = 20;
    int i, j, x, y;
    int n = matrixSize;
    #pragma omp parallel for private(j, x, y)
    for (i = 0; i < n; i += BLOCK_SIZE)
    {
        for (int nn = 0; nn < BLOCK_SIZE; nn++)
        {
            results[i + nn] = 0;
        }

        for (j = 0; j < n; j += BLOCK_SIZE)
        {

            for (x = i; x < fmin(i + BLOCK_SIZE, n); x++)
            {
                #pragma omp simd
                for (y = j; y < (int)fmin(j + BLOCK_SIZE, n); y++)
                {
                    results[x] += A[x*matrixSize + y] * V[y];
                }
            }
        }
    }
    
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

void genRandMatrix(double *A, unsigned long size)
{
    srand(time(0));
    unsigned long i, j;
        for (i = 0; i < size; i++)
        {
            for (j = 0; j < size; j++)
            {
                double n = rand() % 3;
                A[i*size + j] = n;
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

    for (int m = 0; m < NUM_MATRIX_SIZES; m++){

        unsigned long matrixSize = matrixSizes[m]; //

        double *V = (double *)malloc(matrixSize * sizeof(double));
        double *seqV = (double *)malloc(matrixSize * sizeof(double)); // Sequentially computed vector
        double *parV = (double *)malloc(matrixSize * sizeof(double)); // Parallel computed vector

        double *A = (double *)malloc(matrixSize * matrixSize * sizeof(double)); // Matrix to multiply by V



        genRandVector(V, matrixSize);
        genRandMatrix(A, matrixSize);

        for (int a = 0; a < NUM_AVERAGES; a++) {
            clock_gettime(CLOCK_MONOTONIC, &start);
            doSequentialComputation(A, V, seqV, matrixSize);
            clock_gettime(CLOCK_MONOTONIC, &finish);
            elapsed = (finish.tv_sec - start.tv_sec);
            elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
            sequentialTimings[m] += elapsed;
            for (int t = 1; t <= THREAD_RANGE; t++)
            {
                clock_gettime(CLOCK_MONOTONIC, &start);

                doParallelComputation(A, V, parV, matrixSize, t);
                clock_gettime(CLOCK_MONOTONIC, &finish);
                elapsed = (finish.tv_sec - start.tv_sec);
                elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
                parallelTimings[m][t - 1] += elapsed;
                for (int i = 0; i < matrixSize; i++)
                {
                    assert(seqV[i] == parV[i]);
                }
            }
        }
        free(seqV);
        free(parV);
        free(A);
        free(V);
    }

    for (int m = 0; m < NUM_MATRIX_SIZES; m++){
        sequentialTimings[m] /= NUM_AVERAGES;
        for (int t = 0; t < THREAD_RANGE; t++) {
            parallelTimings[m][t] /= NUM_AVERAGES;
        }
    }

    printf("Sequential:\n");
    for (int m = 0; m < NUM_MATRIX_SIZES; m++)
    {
        printf("%lux%lu:\t", matrixSizes[m], matrixSizes[m]);

        printf("%.9f \n",sequentialTimings[m]);
    }

    printf("\nParallel:\n");
    for (int m = 0; m < NUM_MATRIX_SIZES; m++)
    {
        printf("%lux%lu:\t", matrixSizes[m], matrixSizes[m]);

        for (int t = 0; t < THREAD_RANGE; t++)
        {
            printf(" %.9f ",parallelTimings[m][t]);
        }
        printf("\n");
        
    }
    printf("\n");


    return 0;
}