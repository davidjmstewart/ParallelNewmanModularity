#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
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
                int n = rand() % 100;
                if (n < PROBABILITY_OF_CONNECTION){
                    A[i * size + j] = 1;
                    A[j * size + i] = 1;
                } else {
                    A[i * size + j] = 0;
                    A[j * size + i] = 0;
                }

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
        double n = rand() % 3;
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


void printDoubleMatrix(double *matrix, int size)
{
    printf("Doing\n");
    int i, j;
    for (i = 0; i < size; ++i)
    {
        for (j = 0; j < size; ++j)
            printf("%f ", matrix[i*size + j]);
        printf("\n");
    }
}

double dotProduct(double *restrict A, double *restrict x, unsigned long size)
{
    double sum = 0;
    int i;
    #pragma omp simd reduction(+ \
                            : sum)
    for (i = 0; i < size; i++)
    {
        sum += A[i] * x[i];
    }
    return sum;
}

double rayleighQuotient(double *restrict B, double *restrict eigenVector, unsigned long size, int numThreads){
    double e = 0;
    double *By = (double *)malloc(size * sizeof(double));
    matVectMultiply(B, eigenVector, By, size, numThreads);
    e = dotProduct(By, eigenVector, size) / dotProduct(eigenVector, eigenVector, size);
    free(By);
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
        // fwrite("[\n", sizeof(char), 1, fp);

        for (int i = 0; i < size; i++)
        {
            // fwrite("[", sizeof(char), 1, fp);

            for (int j = 0; j < size; j++)
            {
                fprintf(fp, "%d ", A[i * size + j]);
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

void vectorToFile(double *S, unsigned long size, enum OutputType outputType)
{


    if (outputType == MATLAB)
    {
        FILE *fp;
        fp = fopen("./benchmarking/matlab_testing_files/matlab-vec.txt", "w");

        for (int i = 0; i < size; i++)
        {
            fprintf(fp, "%0.9f", S[i]);
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
}

void communitiesToFile(community_t *communities, enum OutputType outputType)
{


    if (outputType == MATLAB)
    {
        FILE *fp;
        fp = fopen("./benchmarking/matlab_testing_files/matlab-communities.txt", "w");
        community_t *tmp = communities;
       
        while (tmp != NULL)
        {
            for (int i = 0; i < tmp->numNodes; i++)
            {
                fprintf(fp, "%d ", tmp->globalVertices[i]+1);
            }
            fprintf(fp, "\n");
            tmp = tmp->nextCommunity;
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
        }
        results[i] = tmp; // write-only to results, not adding to old value.
    }
}

eigenPair powerIteration(double *restrict B, unsigned long size, int numThreads, double tolerance, int iterationLimit)
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

    bool converged = false;
    int numIterations = 0;
    while (!converged && numIterations < iterationLimit)
    {
        int k;
        matVectMultiply(B, eigenVectorTmp, eigenVector, size, numThreads);
        norm = 0;
        #pragma omp simd reduction(+:norm)
        for (k = 0; k < size; k++)
        {
             norm += eigenVector[k] * eigenVector[k];
        }
        // #pragma omp parallel for private(k)
        for (k = 0; k < size; k++)
        {
            eigenVectorTmp[k] = eigenVector[k] / sqrt(norm);
        }
        
        eigenValue = rayleighQuotient(B, eigenVector, size, numThreads);
        double diff = fabs(eigenValue - prevEigenValue);
        if (diff < tolerance)
        {
            // printf("converged after %d iterations, matrix size: %ld \t tolerance: %0.9f\n", numIterations, size, tolerance);
            converged = true;
        }

        prevEigenValue = eigenValue;
        numIterations++;
        
        // numIterations++;
    }

    // N.B. This is not part of the typical power iteration algorithm
    // this is specifically for our purposes of finding the eigenvector
    // corresponding to the most positive eigenvalue. If the eigenvalue is negative we need
    // to perform a spectral shift and repeat the process one more time
    // see these threads:
    //     * https://math.stackexchange.com/questions/835450/efficient-method-for-determining-to-the-most-positive-eigenvalue-of-a-matrix
    //     * https://math.stackexchange.com/questions/906563/finding-eigenvectors-for-the-largest-eigenvalue-vs-one-with-the-largest-absolute

    // double eigenValue = rayleighQuotient(B, eigenVector, size, numThreads);
    if (eigenValue < 0) {
        double *newB = (double *)malloc(size * size * sizeof(double));

        memcpy(newB, B, size * size * sizeof(double));
        for (int i = 0; i < size; i++){
            newB[i * size + i] += fabs(eigenValue); 
        }
        free(eigenVectorTmp);
        free(eigenVector);

        eigP = powerIteration(newB, size, numThreads, tolerance, iterationLimit);

        free(newB);
        return eigP;
        
    } else {
        eigP.eigenvalue = eigenValue;
        eigP.eigenvector = eigenVector;
        return eigP;
    }
    eigP.eigenvalue = eigenValue;
    eigP.eigenvector = eigenVector;
    free(eigenVectorTmp);
    return eigP;
}

void createGlobalVertices(int *restrict globalVertices, unsigned long matrixSize, int numThreads)
{
    for (int i = 0; i < matrixSize; i++)
    {
        globalVertices[i] = i;
    }
}

community_t* assignCommunity(double *restrict B, unsigned long currentMatrixSize, unsigned long originalMatrixSize, int *globalVertices, int numThreads)
{

    double *B_g = (double *)malloc(currentMatrixSize * currentMatrixSize * sizeof(double)); // enough space for a size x size matrix of double types
    struct timespec start, finish;

    computeSubgraphModularityMatrix(B_g, B, currentMatrixSize, originalMatrixSize, globalVertices, numThreads);
    // printf("Doing power iteration in assignCommunity for size %ld \n", currentMatrixSize);
    clock_gettime(CLOCK_MONOTONIC, &start);


    eigenPair eigP = powerIteration(B_g, currentMatrixSize, numThreads, 0.0001, 50);
    clock_gettime(CLOCK_MONOTONIC, &finish);
    double elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    // printf("Took %0.9f seconds for  %ld elements\n", elapsed, currentMatrixSize);

    double *S = (double *)malloc(currentMatrixSize * sizeof(double)); // Membership vector, referred to as S in the paper

    createMembershipVector(eigP.eigenvector, S, currentMatrixSize);

    double *B_gTimesS = (double *)malloc(currentMatrixSize * sizeof(double)); // Store the matrix-vector multiplication of B_g * S here

    matVectMultiply(B_g, S, B_gTimesS, currentMatrixSize, numThreads);
    double deltaQ = 0;
    deltaQ = dotProduct(S, B_gTimesS, currentMatrixSize);
    // printf("DeltaQ: %lf \n", deltaQ);
    if (deltaQ < 0.1 || currentMatrixSize <= 1){ // if this isn't a good split
        // printf("Bad split!\n");
        community_t * community = (community_t *) malloc(sizeof(community_t));
        community->nextCommunity  = NULL;
        community->globalVertices = globalVertices;
        community->numNodes = currentMatrixSize;
        free(B_gTimesS);
        free(B_g);
        free(S);
        // printf("\n");
        return community;
    } else { // if this is a good split

        int *leftIndices =  (int *)malloc(currentMatrixSize * sizeof(int)); // Indices that correspond to values of -1 in S
        int *rightIndices = (int *)malloc(currentMatrixSize * sizeof(int)); // Indices that correspond to values of 1 in S
        int *globalVerticesLeft = (int *)malloc(currentMatrixSize * sizeof(int));  // Global vertices that correspond to values of -1 in S
        int *globalVerticesRight = (int *)malloc(currentMatrixSize * sizeof(int)); // Global vertices that correspond to values of 1 in S
        int numLeft = 0;
        int numRight = 0;

        for (int i = 0; i < currentMatrixSize; i++)
        {
            if (S[i] == 1) {
                rightIndices[numRight] = i;
                globalVerticesRight[numRight] = globalVertices[i];
                numRight++;
            }
            else {
                leftIndices[numLeft] = globalVertices[i];
                globalVerticesLeft[numLeft] = globalVertices[i];
                numLeft++;
            }
        }


        community_t *leftCommunity = malloc(sizeof(struct community));
        community_t *rightCommunity = malloc(sizeof(struct community));

        leftCommunity = assignCommunity(B, numLeft, originalMatrixSize, globalVerticesLeft, numThreads);
        leftCommunity->deltaQ = deltaQ;

        rightCommunity = assignCommunity(B, numRight, originalMatrixSize, globalVerticesRight, numThreads);
        rightCommunity->deltaQ = deltaQ;

        community_t* tmp = leftCommunity;
        while(tmp->nextCommunity != NULL)
        {
            tmp = tmp->nextCommunity;                 // Move to next node
        }
        tmp->nextCommunity = rightCommunity;
        free(leftIndices);
        free(rightIndices);
        free(B_gTimesS);
        free(B_g);

        // free(leftB);
        free(S);

        // free(rightB);
        return leftCommunity;


    }

}

void computeSubgraphModularityMatrix(double *restrict B_g, double *restrict B, unsigned long currentMatrixSize, unsigned long originalMatrixSize, int *globalVertices, int numThreads)
{
    omp_set_num_threads(numThreads);

    // memcpy(B_g, B, currentMatrixSize * currentMatrixSize * sizeof(double));

    for (int i = 0; i < currentMatrixSize; i++)
    {
        for (int j = 0; j < currentMatrixSize; j++)
        {
            B_g[i * currentMatrixSize + j] = B[globalVertices[i] * originalMatrixSize + globalVertices[j]];
        }
    }


    double *Bsum = (double *)malloc(currentMatrixSize * sizeof(double)); // Sum each row of B and put it in here
    int i;

// double tmp;
#pragma omp parallel for private(i)
    for (i = 0; i < currentMatrixSize; i++)
    {
        double *Bhead;
        Bhead = &B_g[i * currentMatrixSize];
        double tmp = 0;
    #pragma omp simd reduction(+ \
                            : tmp)
        for (int j = 0; j < currentMatrixSize; j++)
        {
            tmp += Bhead[j];
        }
        Bsum[i] = tmp;
        // printf("%.4f\n", Bsum[i]);
    }

// #pragma omp parallel for
    for (i = 0; i < currentMatrixSize; i++)
    {
        B_g[i * currentMatrixSize + i] -= Bsum[i];
    }
    free(Bsum);
}


void createMembershipVector(double *restrict S, double *restrict newS, unsigned long currentMatrixSize)
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
