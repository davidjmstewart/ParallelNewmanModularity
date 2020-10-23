import numpy as np
from time import time
def power_iteration(A, num_simulations: int):
    # Ideally choose a random vector
    # To decrease the chance that our vector
    # Is orthogonal to the eigenvector
    # b_k = np.array([1, 1, 1, 1, 1])
    b_k = np.random.rand(A.shape[1])
    for _ in range(num_simulations):
        # calculate the matrix-by-vector product Ab
        b_k1 = np.dot(A, b_k)

        # print("b_k1 at step " + str(_))
        # print(b_k1)
        # print("first element before normalizing: " + str(b_k1[0]))

        # calculate the norm
        b_k1_norm = np.linalg.norm(b_k1)
        # print("norm is: " + str(b_k1_norm))
        # re normalize the vector
        b_k = b_k1 / b_k1_norm
        # print("b_k1 renormalised " + str(b_k))
        # print("")
        # eig_val = d


    return b_k

# print(power_iteration(np.array([
#     [0, 1, 0, 0, 0], 
#     [1, 0, 0, 0, 0], 
#     [0, 0, 0, 0, 1], 
#     [0, 0, 0, 0, 1], 
#     [0, 0, 1, 1, 0], 
#     ]), 10))

B = np.loadtxt('./py-mat.txt')
start = time()
eigen_vector = power_iteration(B, 500)
end = time()

print(end-start)
print("EIGENVECTOR: ")
# print(eigen_vector)
print("By: ")

By = np.matmul(B, eigen_vector)
# print(By)
eigen_value = np.dot(B, eigen_vector) / np.dot(eigen_vector, eigen_vector)
# print(np.dot(np.matmul(B, eigen_vector), eigen_vector))