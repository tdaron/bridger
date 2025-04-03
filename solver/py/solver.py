import numpy as np
from numpy.typing import NDArray

import scipy

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


def c_grad_full(A: NDArray[np.float64], B: NDArray[np.float64], X: NDArray[np.float64], epsilon: float) -> NDArray[np.float64]:
    r = B - A @ X
    d = r.copy()
    r_norm = np.linalg.norm(r)
    if r_norm < epsilon:
        return X
    num_iter = 0
    while True:
        Ad = A @ d
        alpha = r_norm / (d @ Ad)
        X += alpha * d
        if num_iter % 50 == 0:
            r = B - A @ X
        else:
            r -= alpha * Ad
        r_norm_new = r @ r
        if r_norm_new < epsilon:
            return X
        beta = r_norm_new / r_norm
        d = r + beta * d
        r_norm = r_norm_new
        num_iter += 1
        print(f"Iteration {num_iter}, residual norm: {r_norm}")

def pre_c_grad_full(A: NDArray[np.float64], B: NDArray[np.float64], X: NDArray[np.float64], epsilon: float) -> NDArray[np.float64]:
    L, U = ILU(A)
    # C = incomplete_cholesky(A)
    x = np.copy(X)
    r = B - A @ x
    r_norm = np.linalg.norm(r)
    if r_norm < epsilon:
        return x
    z = np.linalg.solve(L, r)
    z = np.linalg.solve(U, z)
    # z = np.linalg.solve(C, r)
    # z = np.linalg.solve(C.T, z)
    d = z.copy()
    rz = r @ z
    num_iter = 0
    while True:
        Ad = A @ d
        alpha = rz / (d @ Ad)
        x += alpha * d
        # Optional recomputation of residual every 50 iterations
        # if num_iter % 50 == 0:
        #     r = B - A @ x
        # else:
        #     r -= alpha * Ad
        r -= alpha * Ad
        r_norm = np.linalg.norm(r)
        if r_norm < epsilon:
            return x
        z = np.linalg.solve(L, r)
        z = np.linalg.solve(U, z)
        # z = np.linalg.solve(C, r)
        # z = np.linalg.solve(C.T, z)
        new_rz = r @ z
        beta = new_rz / rz
        d = z + beta * d
        num_iter += 1
        rz = new_rz
        print(f"Iteration {num_iter}, residual norm: {r_norm}")

def ILU(A: NDArray[np.float64], drop_tol: float = 1e-12) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    n = A.shape[0]
    L = np.identity(n, dtype=np.float64)
    U = np.copy(A)

    for k in range(n):
        for i in range(k+1, n):
            if abs(A[i, k]) > drop_tol:
                L[i, k] = U[i, k] / U[k, k]
                U[i, k] = 0  # Set the lower part of U to zero

                # Only update non-zero elements in the original matrix
                for j in range(k+1, n):
                    if abs(A[i, j]) > drop_tol:
                        U[i, j] = U[i, j] - L[i, k] * U[k, j]
    return L, U

def cholesky(A: NDArray[np.float64]):
    # C = np.zeros_like(A)
    # Extract lower triangular part of A
    C = np.tril(A)
    for i in range(len(C)):
        for j in range(i):
            sum = 0.0
            for k in range(j):
                sum += C[i, k] * C[j, k]
            C[i, j] = (C[i, j] - sum) / C[j, j]
        sum = 0.0
        for k in range(i):
            sum += C[i, k] * C[i, k]
        C[i, i] = np.sqrt(C[i, i] - sum)
    return C

def jacobi_preconditioner(A: NDArray[np.float64]):
    """Create a Jacobi (diagonal) preconditioner for matrix A."""
    diag = np.diag(A).copy()

    # Ensure diagonal elements are non-zero
    diag[np.abs(diag) < 1e-14] = 1.0

    def solve(r: NDArray[np.float64]) -> NDArray[np.float64]:
        return r / diag

    return solve

def incomplete_cholesky(A: NDArray[np.float64]):
    # C = np.zeros_like(A)
    # Extract lower triangular part of A
    C = np.tril(A)
    for i in range(len(A)):
        for j in range(i):
            if A[i, j] == 0.0:
                continue
            sum = 0.0
            for k in range(j):
                sum += C[i, k] * C[j, k]
            C[i, j] = (C[i, j] - sum) / C[j, j]
        if A[i, i] == 0.0:
            continue
        sum = 0.0
        for k in range(i):
            sum += C[i, k] * C[i, k]
        C[i, i] = np.sqrt(C[i, i] - sum)
    return C
