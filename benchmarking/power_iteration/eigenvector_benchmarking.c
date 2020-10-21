#include <time.h>
#include <stdio.h>
#include <math.h>

#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <omp.h>
#include <assert.h>
#include <string.h>
#define N 40000 // Matrix size will be N x N
#define T 1
#define THREAD_RANGE 8 // Run for 1:THREAD_RANGE threads
#define NUM_AVERAGES 20 // take the average of 5 timings for each matrix size, and each number of threads
#define NUM_MATRIX_SIZES 6
#define ITERATION_LIMIT 6000
#include "../../lib/CDUtils.h";
// unsigned long matrixSizes[NUM_MATRIX_SIZES] = {51200};

// unsigned long matrixSizes[NUM_MATRIX_SIZES] = {64, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200};
// unsigned long matrixSizes[NUM_MATRIX_SIZES] = {100000, 200000, 500000, 1000000, 2000000, 3000000, 5000000, 10000000, 20000000, 50000000};
unsigned long matrixSizes[NUM_MATRIX_SIZES] = {128, 256, 512, 1024, 2048, 4096, 8192, 16384};
// unsigned long matrixSizes[NUM_MATRIX_SIZES] = { 51200 };

double sequentialTimings[NUM_MATRIX_SIZES];
double parallelTimings[NUM_MATRIX_SIZES][THREAD_RANGE];
const double TOLERANCE = 0.000000001;
//gcc -fopenmp -g ./benchmarking/matrix_vector_multiplication/mvm_benchmarking.c -o ./benchmarking/matrix_vector_multiplication/mvm ./lib/CDUtils.o -lm -Wall -Wpedantic -Waggressive-loop-optimizations -O3 -march=native -fopt-info-vec-missed=v.txt
// #pragma GCC option("arch=native", "tune=native", "no-zero-upper") //Enable AVX
#pragma GCC target("avx")                                         //Enable AVX

// todo: parallelise
double seqDotProduct(double *A, double *x, unsigned long size)
{
    // double *y = (double *)malloc(size * sizeof(double));
    double sum = 0;
    int i;

    for (i = 0; i < size; i++)
    {
        sum += A[i] * x[i];
    }
    return sum;
}



void seqMatVectMultiply(double *restrict A, double *restrict V, double *restrict results, unsigned long matrixSize)
{
    unsigned long i, j;

    for (i = 0; i < matrixSize; i++)
    {
        double tmp = 0;
        for (j = 0; j < matrixSize; j++)
        {
            tmp += A[i * matrixSize + j] * V[j];
        }
        results[i] = tmp; // write-only to results, not adding to old value.
    }
}

double seqRayleighQuotient(double *B, double *eigenVector, unsigned long size)
{
    double e = 0;
    double *By = (double *)malloc(size * sizeof(double));
    seqMatVectMultiply(B, eigenVector, By, size);
    e = seqDotProduct(By, eigenVector, size) / seqDotProduct(eigenVector, eigenVector, size);
    free(By);
    return e;
}



eigenPair seqPowerIteration(double *B, unsigned long size, double tolerance, int iterationLimit)
{
    double *eigenVectorTmp = (double *)malloc(size * sizeof(double));
    double *eigenVector = (double *)malloc(size * sizeof(double));
    eigenPair eigP;
    genRandMembershipVector(eigenVectorTmp, size);
    for (int i = 0; i < size; i++)
    {
        eigenVector[i] = eigenVectorTmp[i];
    }

    double eigenValue = 0;
    double prevEigenValue = __INT_MAX__;
    double norm = 0;

    bool converged = false;
    int numIterations = 0;
    while (!converged && numIterations < iterationLimit)
    // while (numIterations < iterationLimit)
    {

        int k;

        seqMatVectMultiply(B, eigenVectorTmp, eigenVector, size);

        norm = 0;
// double max = (eigenVector[0]);

        for (k = 0; k < size; k++)
        {
            norm += eigenVector[k] * eigenVector[k];
        }

        for (k = 0; k < size; k++)
        {
            eigenVectorTmp[k] = eigenVector[k] / sqrt(norm);
        }
        eigenValue = seqRayleighQuotient(B, eigenVector, size);
        
        // eigenValue = rayleighQuotient(B, eigenVector, size, numThreads);
        double diff = fabs(eigenValue - prevEigenValue);
        if (diff < tolerance)
        {
            // printf("converged after %d iterations, matrix size: %d\n", numIterations, size);
            converged = true;
        }

        prevEigenValue = eigenValue;
        numIterations++;
        
        // numIterations++;
    }
    // double eigenValue = rayleighQuotient(B, eigenVector, size);

    // N.B. This is not part of the typical power iteration algorithm
    // this is specifically for our purposes of finding the eigenvector
    // corresponding to the most positive eigenvalue. If the eigenvalue is negative we need
    // to perform a spectral shift and repeat the process one more time
    // see these threads:
    //     * https://math.stackexchange.com/questions/835450/efficient-method-for-determining-to-the-most-positive-eigenvalue-of-a-matrix
    //     * https://math.stackexchange.com/questions/906563/finding-eigenvectors-for-the-largest-eigenvalue-vs-one-with-the-largest-absolute

    // double eigenValue = seqRayleighQuotient(B, eigenVector, size);
    
    if (eigenValue < 0)
    {

        double *newB = (double *)malloc(size * size * sizeof(double));

        memcpy(newB, B, size * size * sizeof(double));
        for (int i = 0; i < size; i++)
        {
            newB[i * size + i] += fabs(eigenValue); // todo: change this, it is mutating our original B
        }
        free(eigenVectorTmp);
        free(eigenVector);

        eigP = seqPowerIteration(newB, size, tolerance, iterationLimit);

        free(newB);
        return eigP;
    }
    else
    {
        // printf("eigenvalue is %f\n", eigenValue);
        eigP.eigenvalue = eigenValue;
        eigP.eigenvector = eigenVector;
        return eigP;
        // return eigenVector;
    }
    eigP.eigenvalue = eigenValue;
    eigP.eigenvector = eigenVector;
    free(eigenVectorTmp);
    return eigP;
    // return eigenVector;
}

eigenPair doSequentialComputation(double *B, unsigned long size, double tolerance, int iterationLimit)
{
    return seqPowerIteration(B, size, tolerance, iterationLimit);

}

eigenPair doParallelComputation(double *B, unsigned long size, int numThreads, double tolerance, int iterationLimit)
{

    return powerIteration(B, size, numThreads, tolerance, iterationLimit);




/*
    omp_set_num_threads(numThreads);
    const int BLOCK_SIZE = 256;
    int i, j, x, y;
    int n = matrixSize;
#pragma omp parallel for private(j, x, y)
    for (i = 0; i < n; i += BLOCK_SIZE)
    {
        for (int nn = 0; nn < BLOCK_SIZE; nn++)
        {
            results[i + nn] = 0;
        }
        int xmin = (i + BLOCK_SIZE < n ? i + BLOCK_SIZE : n);
        for (j = 0; j < n; j += BLOCK_SIZE)
        {
            int ymin = (j + BLOCK_SIZE < n ? j + BLOCK_SIZE : n);
            for (x = i; x < xmin; x++)
            {
                double tmp = 0;
                double *MHead = &M[x * matrixSize];
#pragma omp simd reduction(+ \
                           : tmp)
                for (y = j; y < ymin; y++)
                {
                    tmp += MHead[y] * V[y];
                }
                results[x] += tmp;
            }
        }
    }
    
*/
    
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
            eigenPair eigPairSequential = doSequentialComputation(A, matrixSize, TOLERANCE, ITERATION_LIMIT);
            clock_gettime(CLOCK_MONOTONIC, &finish);
            elapsed = (finish.tv_sec - start.tv_sec);
            elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
            sequentialTimings[m] += elapsed;
            for (int t = 1; t <= THREAD_RANGE; t++)
            {
                clock_gettime(CLOCK_MONOTONIC, &start);

                eigenPair eigPairParallel = doParallelComputation(A, matrixSize, t, TOLERANCE, ITERATION_LIMIT);
                // printf("Seq: %f \t Par: %f \n", eigPairSequential.eigenvalue, eigPairParallel.eigenvalue);
                clock_gettime(CLOCK_MONOTONIC, &finish);
                elapsed = (finish.tv_sec - start.tv_sec);
                elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
                parallelTimings[m][t - 1] += elapsed;
                assert(fabs(eigPairParallel.eigenvalue - eigPairSequential.eigenvalue) < 0.1);
                // for (int i = 0; i < matrixSize; i++)
                // {
                //     assert(seqV[i] == parV[i]);
                // }
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