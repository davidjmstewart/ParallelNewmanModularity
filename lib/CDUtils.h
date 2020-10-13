#ifndef CDUTILS_H
#define CDUTILS_H

#define USE_DEFAULT_NUM_THREADS -1
#define DEFAULT_NUM_THREADS 16

enum OutputType
{
    Python = 1,
    MATLAB = 0
};

typedef struct eigenPair
{
    double eigenvalue;
    double *eigenvector;
} eigenPair;

// node for a linked list
typedef struct node
{
    int num;
    struct node *next;
} node_t;

// Generates a size x size adjacency matrix in a flat structure
// INDEXING EXAMPLE
// If size is 5, the 3rd row, 4th column of A would be A[2*size + 3]
void genAdjacencyMatrix(int *A, unsigned long size);

// Creates a vector that stores the degree of each node in A
// Takes a (size x size) adjacency matrix A (e.g. one made with genAdjacencyMatrix())
// and places the sum of each row into vector D of length size
void createDegreesVec(int *A, int *D, unsigned long size, int numThreads);

// Computes equation 3 in the paper. B is the modularity matrix that will be filled
// by this function, A is the symmetric adjacency matrix describing this graph
// D is the degree vector where D[i] contains the degree of A[i, :]
void createModularityMatrix(double *B, int *A, int *D, unsigned long size, int numThreads);

int graphDegree(int *D, unsigned long size);

eigenPair powerIteration(double *B, unsigned long size, int numThreads, double tolerance, int iterationLimit);

void matrixToFile(double *A, unsigned long size, enum OutputType outputType);
void membershipVectorToFile(double *S, unsigned long size, enum OutputType outputType);
void integerMatrixToFile(int *A, unsigned long size, enum OutputType outputType);

void matVectMultiply(double *restrict A, double *restrict V, double *restrict results, unsigned long matrixSize, int numThreads);

// bulk of Newman's algorithm driver logic is in here.
void assignCommunity(double *restrict B, int nextGroupNum, unsigned long currentMatrixSize, unsigned long originalMatrixSize, int *globalVertices, int numThreads);

// Computes equation 6 in the paper. Places the results in B_g
void computeSubgraphModularityMatrix(double *restrict B_g, double *restrict B, unsigned long currentMatrixSize, unsigned long originalMatrixSize, int *globalVertices, int numThreads);

#endif /* CDUTILS_H */