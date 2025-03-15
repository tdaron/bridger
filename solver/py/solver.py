import numpy as np
from numpy.typing import NDArray

def np_solve(A_global: NDArray[np.float64], B_global: NDArray[np.float64]) -> NDArray[np.float64]:
    U = np.linalg.solve(A_global, B_global)
    return U.astype(np.float64)
