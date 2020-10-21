import os


# os.environ['OPENBLAS_NUM_THREADS'] = '1'
# os.environ['MKL_NUM_THREADS'] = '1'
# os.environ["MKL_NUM_THREADS"] = "1" 
# os.environ["NUMEXPR_NUM_THREADS"] = "1" 
# os.environ["OMP_NUM_THREADS"] = "1" 
import numpy as np
from numpy import linalg as LA
import sys
import time
import timeit 

NUM_AVERAGES = 10

average_timings = []

matrix_sizes = [128, 256, 512, 1024, 2048, 4096, 8192, 16384]
def time_computation(size, num_averages):
    A = np.random.rand(size, size)
    V = np.random.rand(size, 1)

    def closure():
        return np.multiply(A, V)
    return timeit.timeit(closure, number=num_averages)

print(time_computation(100, 1000))


for i, size in enumerate(matrix_sizes):
    A = np.random.rand(size, size).astype('f')
    V = np.random.rand(size, 1).astype('f')
    size_avg = 0
    for avg_iter in range(NUM_AVERAGES):
        # start_time = time.monotonic()
        size_avg = time_computation(size, NUM_AVERAGES)/NUM_AVERAGES


    average_timings.append(size_avg)

print(average_timings)