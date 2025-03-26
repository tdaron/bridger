import numpy as np
from scipy.sparse import csr_matrix
import time

class Problem:
    """Mock Problem class for testing purposes"""
    def __init__(self, elements, points, n_local):
        self.elements = elements
        self.points = points
        self.n_local = n_local

def rcm(problem: Problem) -> tuple:
    """Your implementation of CSR matrix creation"""
    # Initialize CSR matrix lists
    coo_adj = []
    for i in range(problem.elements.shape[0]):
        for j in range(problem.n_local):
            for k in range(problem.n_local):
                coo_adj.append((problem.elements[i, j], problem.elements[i, k]))

    # Sort the list of tuples
    coo_adj.sort(key=lambda x: (x[0], x[1]))

    # Create the CSR matrix
    col_index = []
    row_index = list(np.zeros(problem.points.shape[0] + 1, dtype=int))
    row_index = [0 for _ in range(problem.points.shape[0] + 1)]

    prev_row = coo_adj[0][0]
    prev = None
    for i in range(len(coo_adj)):
        # Skip duplicates
        if prev == coo_adj[i]:
            continue
        col_index.append(coo_adj[i][1])
        # Set all absent rows that were skipped to the
        # current column index.
        for j in range(prev_row, coo_adj[i][0]):
            row_index[j + 1] = len(col_index) - 1
        prev_row = coo_adj[i][0]
        prev = coo_adj[i]

    # Make sure to update the last rows that might be empty
    for j in range(prev_row + 1, len(row_index)):
        row_index[j] = len(col_index)

    return col_index, row_index

def scipy_csr_matrix(problem: Problem) -> csr_matrix:
    """Create CSR matrix using SciPy for comparison"""
    rows = []
    cols = []
    data = []

    for i in range(problem.elements.shape[0]):
        for j in range(problem.n_local):
            for k in range(problem.n_local):
                rows.append(problem.elements[i, j])
                cols.append(problem.elements[i, k])
                data.append(1.0)  # Using 1.0 as placeholder value

    return csr_matrix((data, (rows, cols)), shape=(problem.points.shape[0], problem.points.shape[0]))

def verify_csr_equivalence(problem: Problem) -> bool:
    """Verify if your CSR implementation is equivalent to SciPy's"""
    # Get your implementation result
    col_index, row_index = rcm(problem)

    # Get SciPy's implementation
    scipy_matrix = scipy_csr_matrix(problem)

    # Extract SciPy's CSR components
    scipy_indices = scipy_matrix.indices
    scipy_indptr = scipy_matrix.indptr

    # Check if the structures match
    indices_match = np.array_equal(col_index, scipy_indices)
    indptr_match = np.array_equal(row_index, scipy_indptr)

    if not indices_match:
        print("Column indices don't match:")
        print("Your implementation:", col_index)
        print("SciPy implementation:", scipy_indices)

    if not indptr_match:
        print("Row pointers don't match:")
        print("Your implementation:", row_index)
        print("SciPy implementation:", scipy_indptr)

    return indices_match and indptr_match

def test_rcm():
    """Test the RCM implementation with different cases"""

    # Test Case 1: Simple matrix with no gaps
    elements = np.array([[0, 1], [1, 2], [2, 0]])
    points = np.zeros((3, 2))  # 3 points
    problem1 = Problem(elements, points, 2)

    # Test Case 2: Matrix with gaps (missing rows)
    elements = np.array([[0, 3], [3, 6], [6, 9]])
    points = np.zeros((10, 2))  # 10 points, but only indices 0,3,6,9 are used
    problem2 = Problem(elements, points, 2)

    # Test Case 3: Larger random matrix
    n_elem = 100
    n_points = 200
    n_local = 3
    elements = np.random.randint(0, n_points, (n_elem, n_local))
    points = np.zeros((n_points, 2))
    problem3 = Problem(elements, points, n_local)

    # Run tests
    print("Test Case 1:", "PASSED" if verify_csr_equivalence(problem1) else "FAILED")
    print("Test Case 2:", "PASSED" if verify_csr_equivalence(problem2) else "FAILED")
    print("Test Case 3:", "PASSED" if verify_csr_equivalence(problem3) else "FAILED")

    # Performance test
    start = time.time()
    rcm(problem3)
    end = time.time()
    print(f"Your implementation time: {end - start:.6f} seconds")

    start = time.time()
    scipy_csr_matrix(problem3)
    end = time.time()
    print(f"SciPy implementation time: {end - start:.6f} seconds")

if __name__ == "__main__":
    test_rcm()
