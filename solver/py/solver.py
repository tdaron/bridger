import numpy as np
from numpy.typing import NDArray

def np_solve(A_global: NDArray[np.float64], B_global: NDArray[np.float64]) -> NDArray[np.float64]:
    U = np.linalg.solve(A_global, B_global)
    return U.astype(np.float64)

def sym_band_cholesky(A: NDArray[np.float64], B: NDArray[np.float64], U: NDArray[np.float64], bandwidth: int):
    """
    Symmetric Gaussian elimination for banded matrices
    """
    idx = lambda i, j: (i + 1) * bandwidth + j

    size = len(B)

    # Cholesky decomposition
    for i in range(size):
        j0 = max(0, i - bandwidth)
        for j in range(j0, i):
            sum = 0.0
            for k in range(j0, j):
                sum += A[idx(i, k)] * A[idx(j, k)]
            A[idx(i, j)] = (A[idx(i, j)] - sum) / A[idx(j, j)]
        sum = 0.0
        for k in range(j0, i):
            sum += A[idx(i, k)] * A[idx(i, k)]
        A[idx(i, i)] = np.sqrt(A[idx(i, i)] - sum)

    # Forward substitution
    # Cy = b
    for i in range(size):
        j0 = max(0, i - bandwidth)
        sum = 0.0
        for j in range(j0, i):
            sum += A[idx(i, j)] * U[j]
        U[i] = (B[i] - sum) / A[idx(i, i)]

    # Backward substitution
    # C^Tx = y
    for i in range(size - 1, -1, -1):
        j0 = min(size, i + bandwidth + 1)
        sum = 0.0
        for j in range(i + 1, j0):
            sum += A[idx(j, i)] * U[j]
        U[i] = (U[i] - sum) / A[idx(i, i)]

    return U.astype(np.float64)
