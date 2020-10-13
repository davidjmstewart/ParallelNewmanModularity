#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

#include "./CDUTils.h"

// gcc -o CDUtils.o -c CDUtils.c
// Output assembly: gcc -fopenmp -O3 -march=native -S CDUtils.c
// Generates a size x size adjacency matrix in a flat structure
// INDEXING EXAMPLE
// If size is 5, the 3rd row, 4th column of A would be A[2*size + 3]
void genAdjacencyMatrix(int *A, unsigned long size) {
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
    
    #pragma omp parallel for private(i, j) 
    for (i = 0; i < size; i++)
    {
        int sum = 0;
        int *ARow = &A[i * size]; // compute row once in this thread to save constant dereferencing
        for (j = 0; j < size; j++)
        {
            sum += *(ARow + j);
        }
        D[i] = sum;
    }
}

// Computes equation 3 in the paper. B is the modularity matrix that will be filled
// by this function, A is the symmetric adjacency matrix describing this graph
// D is the degree vector where D[i] contains the degree of A[i, :]
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

    for (i = 0; i < size; i++)
    {
        sum += A[i] * x[i];
    }
    return sum;
}

double rayleighQuotient(double * B, double * eigenVector, unsigned long size, int numThreads){
    double e = 0;
    double *By = (double *)malloc(size * sizeof(double));
    // int i, j;

    matVectMultiply(B, eigenVector, By, size, numThreads);

    /*
    #pragma omp parallel for private(i, j)
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
    */

    e = dotProduct(By, eigenVector, size) / dotProduct(eigenVector, eigenVector, size);

    return e;
}


void matrixToFile(double *A, unsigned long size, enum OutputType outputType) {
    
    if (outputType == Python)
    {
        FILE *fp = NULL;
        fp = fopen("./benchmarking/python_testing_files/py-mat.txt", "w");
        // fwrite("[\n", sizeof(char), 1, fp);

        for (int i = 0; i < size; i++)
        {
            // fwrite("[", sizeof(char), 1, fp);

            for (int j = 0; j < size; j++)
            {
                fprintf(fp, "%0.9f ", A[i * size + j]);
                // if (j != size - 1)
                // {
                //     fwrite(", ", sizeof(char), 2, fp);
                // }
            }
            fwrite("\n", sizeof(char), 1, fp);
        }
        // fwrite("]\n", sizeof(char), 1, fp);
        fclose(fp);
    }
    else if (outputType == MATLAB)
    {
        FILE *fp = NULL;
        fp = fopen("./benchmarking/matlab_testing_files/matlab-mat.txt", "w");

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
        fclose(fp);
    }

}

void integerMatrixToFile(int *A, unsigned long size, enum OutputType outputType)
{

    if (outputType == Python)
    {
        FILE *fp = NULL;
        fp = fopen("./benchmarking/python_testing_files/integer-matrix-py-mat.txt", "w");
        fwrite("[\n", sizeof(char), 1, fp);

        for (int i = 0; i < size; i++)
        {
            fwrite("[", sizeof(char), 1, fp);

            for (int j = 0; j < size; j++)
            {
                fprintf(fp, "%d", A[i * size + j]);
                if (j != size - 1)
                {
                    fwrite(", ", sizeof(char), 2, fp);
                }
            }
            fwrite("],\n", sizeof(char), 3, fp);
        }
        fwrite("]\n", sizeof(char), 1, fp);
        fclose(fp);
    }
    else if (outputType == MATLAB)
    {
        FILE *fp = NULL;
        fp = fopen("./benchmarking/matlab_testing_files/integer-matrix-matlab-mat.txt", "w");

        for (int i = 0; i < size; i++)
        {

            for (int j = 0; j < size; j++)
            {
                // fwrite(&A[i * size + j], sizeof(double), 1, fp); /* Write to File */
                fprintf(fp, "%d", A[i * size + j]);
                if (j != size - 1)
                {
                    fwrite(" ", sizeof(char), 1, fp);
                }
            }
            fwrite("\n", sizeof(char), 1, fp);
        }
        fclose(fp);
    }
}

void membershipVectorToFile(double *S, unsigned long size, enum OutputType outputType)
{


    if (outputType == Python)
    {
        FILE *fp;
        fp = fopen("./benchmarking/python_testing_files/py-membershipvec.txt", "w");
        for (int i = 0; i < size; i++)
        {
            int output = 1;
            if (S[i] < 0){
                output = -1;
            }
            fprintf(fp, "%d\n", output);
        }
        fclose(fp);
    }
    else if (outputType == MATLAB)
    {
        FILE *fp;
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
        fclose(fp);
    }

}

void matVectMultiply(double *restrict A, double *restrict V, double *restrict results, unsigned long matrixSize, int numThreads)
{
    omp_set_num_threads(numThreads);
    unsigned long i, j;

    #pragma omp parallel for private(j)
    for (i = 0; i < matrixSize; i++)
    {
        double tmp = 0;
    // TODO: unroll outer loop and cache-block it.
    #pragma omp simd reduction(+ \
                           : tmp)
        for (j = 0; j < matrixSize; j ++)
        {
            tmp += A[i * matrixSize + j] * V[j];
            // tmp += A[i * matrixSize + j + 1] * V[j + 1];
            // tmp += A[i * matrixSize + j + 2] * V[j + 2];
            // tmp += A[i * matrixSize + j + 3] * V[j + 3];
        }
        results[i] = tmp; // write-only to results, not adding to old value.
    }
}

eigenPair powerIteration(double *B, unsigned long size, int numThreads, double tolerance, int iterationLimit)
{
    double *eigenVectorTmp = (double *)malloc(size * sizeof(double));
    double *eigenVector = (double *)malloc(size * sizeof(double));
    eigenPair eigP;
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

        /*
        int i, j, k;
        #pragma omp parallel for private(j)
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
    */
        int k;

        matVectMultiply(B, eigenVectorTmp, eigenVector, size, numThreads);

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
        eigenValue = rayleighQuotient(B, eigenVector, size, numThreads);
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
        eigP.eigenvalue = eigenValue;
        eigP.eigenvector = eigenVector;
        return eigP;
        // return eigenVector;
    }
    eigP.eigenvalue = eigenValue;
    eigP.eigenvector = eigenVector;
    return eigP;
    // return eigenVector;
}

void assignCommunity(double *restrict B, int nextGroupNum, unsigned long currentMatrixSize, unsigned long originalMatrixSize, int *globalVertices, int numThreads)
{
    double *B_g = (double *)malloc(currentMatrixSize * currentMatrixSize * sizeof(double)); // enough space for a size x size matrix of double types
    computeSubgraphModularityMatrix(B_g, B, currentMatrixSize, originalMatrixSize, globalVertices, numThreads);
    //
    // double *B_g = (double *)malloc(size * size * sizeof(double)); // enough space for a size x size matrix of double types

    // computeSubgraphModularityMatrix(B_g, B, size, numThreads);
    // double *eigenVector = powerIteration(B_g, size, numThreads, 0.0001, 5000);
    // eigenPair eigP = powerIteration(B_g, size, numThreads, 0.0001, 5000); 
    // double * = (double *)malloc(size * sizeof(double));
    // createMembershipVector(eigenVector, S, unsigned long matrixSize);

    /*
        if (eigenPair.eigenvalue < 0.1){


        } else {

            int numLeft, numRight = 0;

            for (int i = 0; i < size; i++){
                if (S[i] == 1){
                    numRight++;
                } else {
                    numLeft++;
                }   
            }

        }

    */

    /*
#pragma omp task shared(i) firstprivate(n)
    i = fib(n - 1);

#pragma omp task shared(j) firstprivate(n)
    j = fib(n - 2);

#pragma omp taskwait
*/

    // free(S);

    // free(B_g);
}



void computeSubgraphModularityMatrix(double *restrict B_g, double *restrict B, unsigned long currentMatrixSize, unsigned long originalMatrixSize, int *globalVertices, int numThreads)
{
    omp_set_num_threads(numThreads);
#pragma omp parallel for
    for (int i = 0; i < currentMatrixSize; i++)
    {
        for (int j = 0; j < currentMatrixSize; j++)
        {
            B_g[i * currentMatrixSize + j] = B[globalVertices[i] * originalMatrixSize + globalVertices[j]];
        }
    }
    int i;
    // int n = currentMatrixSize;
    double *Bsum = (double *)malloc(currentMatrixSize * sizeof(double)); // Sum each row of B and put it in here

// double tmp;
#pragma omp parallel for
    for (i = 0; i < currentMatrixSize; i++)
    {
        double *Bhead;
        Bhead = &B_g[i * currentMatrixSize];
        double tmp = 0;
        for (int j = 0; j < currentMatrixSize; j++)
        {
            tmp += Bhead[j];
        }
        Bsum[i] = tmp;
        // printf("%.4f\n", Bsum[i]);
    }

#pragma omp parallel for
    for (i = 0; i < currentMatrixSize; i++)
    {
        B_g[i * currentMatrixSize + i] -= Bsum[i];
    }
}

void createMembershipVector(double *restrict S, double *restrict newS, unsigned long currentMatrixSize, unsigned long originalMatrixSize)
{
    for (int i = 0; i < currentMatrixSize; i++)
    {
        if (S[i] < 0)
        {
            newS[i] = 1;
        }
        else
        {
            newS[i] = -1;
        }
    }
}
