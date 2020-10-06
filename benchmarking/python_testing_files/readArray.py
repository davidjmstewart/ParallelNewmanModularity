import numpy as np
from numpy import linalg as LA
import time

data = np.loadtxt('./mat.txt')
t1 = time.time()
w, v = LA.eig(data)
print(w)
print(v)
t2 = time.time()
print(t2-t1)
# print(data)