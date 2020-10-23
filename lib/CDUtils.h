#ifndef CDUTILS_H
#define CDUTILS_H

#define USE_DEFAULT_NUM_THREADS -1
#define DEFAULT_NUM_THREADS 16
#define MIN_MATRIX_TASK_SIZE 3000 // matrix must be at least this size in rows and columns before we use tasks to split up working on it
#define PROBABILITY_OF_CONNECTION 30 // Probability that any two nodes are connected (used when generating adjacency graph)

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

typedef struct community
{
    int numNodes;               // how many nodes are in this community
    int *globalVertices;        // global index of nodes that belong in this community
    double deltaQ;              // the change in modularity score made by this community
    struct community *nextCommunity;
} community_t;

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

void communitiesToFile(community_t *communities, enum OutputType outputType);
void vectorToFile(double *S, unsigned long size, enum OutputType outputType);

// Creates a vector that stores the degree of each node in A
// Takes a (size x size) adjacency matrix A (e.g. one made with genAdjacencyMatrix())
// and places the sum of each row into vector D of length size
void createDegreesVec(int *A, int *D, unsigned long size, int numThreads);

// Computes equation 3 in the paper. B is the modularity matrix that will be filled
// by this function, A is the symmetric adjacency matrix describing this graph
// D is the degree vector where D[i] contains the degree of A[i, :]
void createModularityMatrix(double *B, int *A, int *D, unsigned long size, int numThreads);
void createMembershipVector(double *restrict S, double *restrict newS, unsigned long currentMatrixSize);
int graphDegree(int *D, unsigned long size);

eigenPair powerIteration(double *B, unsigned long size, int numThreads, double tolerance, int iterationLimit);

void matrixToFile(double *A, unsigned long size, enum OutputType outputType);
void membershipVectorToFile(double *S, unsigned long size, enum OutputType outputType);
void integerMatrixToFile(int *A, unsigned long size, enum OutputType outputType);

void matVectMultiply(double *restrict A, double *restrict V, double *restrict results, unsigned long matrixSize, int numThreads);

// bulk of Newman's algorithm driver logic is in here.
community_t* assignCommunity(double *restrict B, unsigned long currentMatrixSize, unsigned long originalMatrixSize, int *globalVertices, int numThreads);
void printDoubleMatrix(double *matrix, int size);
void createGlobalVertices(int *restrict globalVertices, unsigned long matrixSize, int numThreads);

// Computes equation 6 in the paper. Places the results in B_g
void computeSubgraphModularityMatrix(double *restrict B_g, double *restrict B, unsigned long currentMatrixSize, unsigned long originalMatrixSize, int *globalVertices, int numThreads);

#endif /* CDUTILS_H */