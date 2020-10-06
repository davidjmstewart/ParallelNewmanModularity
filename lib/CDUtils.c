#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

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

void genRandMembershipVector(double *S, unsigned long size)
{
    srand(time(0));
    unsigned long i;
    #pragma omp parallel for private(i)
    for ( i = 0; i < size; i++)
    {
        double n = rand() % 2;
        n = (n == 0 ? -1 : 1);
        S[i] = n;
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

// todo: parallelise
double dotProduct(double *A, double *x, unsigned long size)
{
    // double *y = (double *)malloc(size * sizeof(double));
    double sum = 0;
    int i;
    #pragma omp parallel for private(i)
    for (i = 0; i < size; i++)
    {
        sum += A[i] * x[i];
    }
    return sum;
}

double rayleighQuotient(double * B, double * eigenVector, unsigned long size){
    double e = 0;
    double *By = (double *)malloc(size * sizeof(double));
    int i, j;

    for (i = 0; i < size; i++)
    {
        double *BHead = &B[i * size];
        double tmp = 0;

        for (j = 0; j < size; j++)
        {
            tmp += BHead[j] * eigenVector[j];
        }
        By[i] = tmp;
    }


    e = dotProduct(By, eigenVector, size) / dotProduct(eigenVector, eigenVector, size);

    return e;
}


void matrixToFile(double *A, unsigned long size, enum OutputType outputType) {
    FILE *fp;

    if (outputType == Python){
        fp = fopen("./py-mat.txt", "w");
        fwrite("[\n", sizeof(char), 1, fp);

        for (int i = 0; i < size; i++)
        {
            fwrite("[", sizeof(char), 1, fp);

            for (int j = 0; j < size; j++)
            {
                fprintf(fp, "%0.9f", A[i * size + j]);
                if (j != size - 1)
                {
                    fwrite(", ", sizeof(char), 2, fp);
                }
            }
            fwrite("],\n", sizeof(char), 3, fp);
        }
        fwrite("]\n", sizeof(char), 1, fp);
    }
    else if (outputType == MATLAB)
    {
        fp = fopen("./matlab-mat.txt", "w");

        for (int i = 0; i < size; i++)
        {

            for (int j = 0; j < size; j++)
            {
                // fwrite(&A[i * size + j], sizeof(double), 1, fp); /* Write to File */
                fprintf(fp, "%0.9f", A[i * size + j]);
                if (j != size - 1){
                    fwrite(" ", sizeof(char), 1, fp);
                }

            }
            fwrite("\n", sizeof(char), 1, fp);
        }
    }

    fclose(fp);
}

void membershipVectorToFile(double *S, unsigned long size, enum OutputType outputType)
{
    FILE *fp;

    if (outputType == Python)
    {
        fp = fopen("./benchmarking/python_testing_files/py-membershipvec.txt", "w");
        for (int i = 0; i < size; i++)
        {
            int output = 1;
            if (S[i] < 0){
                output = -1;
            }
            fprintf(fp, "%d\n", output);
        }
    }
    else if (outputType == MATLAB)
    {
        fp = fopen("./benchmarking/matlab_testing_files/matlab-membershipvec.txt", "w");

        for (int i = 0; i < size; i++)
        {
            int output = 1;
            if (S[i] < 0)
            {
                output = -1;
            }
            fprintf(fp, "%d\n", output);
        }
    }

    fclose(fp);
}

double* powerIteration(double *B, unsigned long size, int numThreads, double tolerance, int iterationLimit)
{
    double *eigenVectorTmp = (double *)malloc(size * sizeof(double));
    double *eigenVector = (double *)malloc(size * sizeof(double));

    genRandMembershipVector(eigenVectorTmp, size);
    for (int i = 0; i < size; i++){
        eigenVector[i] = eigenVectorTmp[i];
    }

    double eigenValue = 0;
    double prevEigenValue = __INT_MAX__;
    double norm = 0;

    omp_set_num_threads(numThreads);
    bool converged = false;
    int numIterations = 0;
    while (!converged && numIterations < iterationLimit)
    {

        int i, j, k;
        #pragma omp parallel for private(i, j, k)
        for (i = 0; i < size; i++)
        {
            double *BHead = &B[i * size];
            double tmp = 0;
    
            for (j = 0; j < size; j++)
            {
                tmp += BHead[j] * eigenVectorTmp[j];
            }
            eigenVector[i] = tmp;

        }

        norm = 0;
        // double max = (eigenVector[0]);
        #pragma omp parallel for private(k) reduction(+:norm)
        for (k = 0; k < size; k++)
        {
             norm += eigenVector[k] * eigenVector[k];
        }
        #pragma omp parallel for private(k)
        for (k = 0; k < size; k++)
        {
            eigenVectorTmp[k] = eigenVector[k] / sqrt(norm);
        }
        eigenValue = rayleighQuotient(B, eigenVector, size);
        double diff = fabs(eigenValue - prevEigenValue);
        if (diff < tolerance)
        {
            printf("converged after %d iterations \n", numIterations);
            converged = true;
        }

        prevEigenValue = eigenValue;
        numIterations++;
    }
    // double eigenValue = rayleighQuotient(B, eigenVector, size);

    // N.B. This is not part of the typical power iteration algorithm
    // this is specifically for our purposes of finding the eigenvector
    // corresponding to the most positive eigenvalue. If the eigenvalue is negative we need
    // to perform a spectral shift and repeat the process one more time
    // see these threads:
    //     * https://math.stackexchange.com/questions/835450/efficient-method-for-determining-to-the-most-positive-eigenvalue-of-a-matrix
    //     * https://math.stackexchange.com/questions/906563/finding-eigenvectors-for-the-largest-eigenvalue-vs-one-with-the-largest-absolute

    
    if (eigenValue < 0) {
        for (int i = 0; i < size; i++){
            B[i * size + i] += fabs(eigenValue);
        }
        free(eigenVectorTmp);
        free(eigenVector);
        return powerIteration(B, size, numThreads, tolerance, iterationLimit);
    } else {
        printf("eigenvalue is %f\n", eigenValue);
        return eigenVector;
    }
    return eigenVector;
}
