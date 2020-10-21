import os
import numpy as np
from numpy import linalg as LA
import sys
import time
import timeit 

NUM_AVERAGES = 1

average_timings = []
print("Started!")
matrix_sizes = [500]
def time_computation(size, num_averages):
    A = np.random.rand(size, size)

    def closure():
        return  LA.eig(A)

    return timeit.timeit(closure, number=num_averages)


for i, size in enumerate(matrix_sizes):
    size_avg = 0
    for avg_iter in range(NUM_AVERAGES):
        # start_time = time.monotonic()
        size_avg = time_computation(size, NUM_AVERAGES)/NUM_AVERAGES


    average_timings.append(size_avg)

print(average_timings)