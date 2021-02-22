#!/usr/bin/env python3

# modified from https://gist.github.com/markus-beuckelmann/8bc25531b11158431a5b09a45abd6276

import numpy as np
from time import time

# Let's take the randomness out of random numbers (for reproducibility)
np.random.seed(0)

size = 4096
A, B = np.random.random((int(size / 4), int(size / 4))), np.random.random((int(size / 4), int(size / 4)))
C, D = np.random.random((size * 128,)), np.random.random((size * 128,))
E = np.random.random((int(size / 2), int(size / 4)))
F = np.random.random((int(size / 2), int(size / 2)))
F = np.dot(F, F.T)
G = np.random.random((int(size / 2), int(size / 2)))

# Matrix multiplication
N = 100
delta = 0
for i in range(N):
    t = time()
    np.dot(A, B)
    delta += time() - t
print('Dotted two %dx%d matrices in %0.2f ms (N = %d).' % (size / 4, size / 4, 1e3 * delta / N, N))

# Vector multiplication
N = 5000
delta = 0
for i in range(N):
    t = time()
    np.dot(C, D)
    delta += time() - t
print('Dotted two vectors of length %d in %0.2f ms (N = %d).' % (size * 128, 1e3 * delta / N, N))

# Singular Value Decomposition (SVD)
N = 3
delta = 0
for i in range(N):
    t = time()
    np.linalg.svd(E, full_matrices = False)
    delta += time() - t
print("SVD of a %dx%d matrix in %0.2f s (N = %d)." % (size / 2, size / 4, delta / N, N))

# Cholesky Decomposition
N = 10
delta = 0
for i in range(N):
    t = time()
    np.linalg.cholesky(F)
    delta += time() - t
print("Cholesky decomposition of a %dx%d matrix in %0.2f s (N = %d)." % (size / 2, size / 2, delta / N, N))

# Eigendecomposition
N = 1
delta = 0
for i in range(N):
    t = time()
    np.linalg.eig(G)
    delta += time() - t
print("Eigendecomposition of a %dx%d matrix in %0.2f s (N = %d)." % (size / 2, size / 2, delta / N, N))

print('')
print('This was obtained using the following NumPy configuration:')
np.__config__.show()
