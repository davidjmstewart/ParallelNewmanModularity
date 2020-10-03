#ifndef CDUTILS_H
#define CDUTILS_H

#define USE_DEFAULT_NUM_THREADS -1
#define DEFAULT_NUM_THREADS 16

// Generates a size x size adjacency matrix in a flat structure
// INDEXING EXAMPLE
// If size is 5, the 3rd row, 4th column of A would be A[2*size + 3]
void genAdjacenyMatrix(int *A, unsigned long size);

// Creates a vector that stores the degree of each node in A
// Takes a (size x size) adjacency matrix A (e.g. one made with genAdjacencyMatrix())
// and places the sum of each row into vector D of length size
void createDegreesVec(int *A, int *D, unsigned long size, int numThreads);

void createModularityMatrix(double *B, int *A, int *D, unsigned long size, int numThreads);

int graphDegree(int *D, unsigned long size);

#endif /* CDUTILS_H */