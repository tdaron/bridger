from problem import Problem

import numpy as np
from typing import List

# x axis renumbering
def renumber_x_axis(problem: Problem) -> Problem:
    id_coords_list = list(enumerate(problem.points))

    # sort by x coordinate
    id_coords_list.sort(key=lambda x: x[1][0])

    # renumber
    # for edge in problem.edges:
    #     for i in range(2):
    #         edge[i] = [x[0] for x in id_coords_list].index(edge[i])

    # for element in problem.elements:
    #     for i in range(problem.n_local):
    #         element[i] = [x[0] for x in id_coords_list].index(element[i])

    problem.renum = [x[0] for x in id_coords_list]
    problem.order = [x[0] for x in id_coords_list]
    problem.renum = list(np.zeros(len(id_coords_list), dtype=int))
    for i in range(len(id_coords_list)):
        problem.renum[id_coords_list[i][0]] = i
    # problem.order = [x[0] for x in enumerate(id_coords_list)]
    # problem.renum = [x[0] for x in enumerate(id_coords_list)]

    # problem.points = np.array([x[1] for x in id_coords_list])

    return problem

def problem_coo(problem: Problem) -> List[tuple]:
    """Your implementation of CSR matrix creation"""
    # Initialize CSR matrix lists
    coo_adj = []
    for i in range(problem.elements.shape[0]):
        for j in range(problem.n_local):
            for k in range(problem.n_local):
                coo_adj.append((problem.elements[i, j], problem.elements[i, k]))

    # Sort the list of tuples
    coo_adj.sort(key=lambda x: (x[0], x[1]))

    return coo_adj

def csr_from_coo(coo_adj: List[tuple], n_points: int) -> tuple:
    """Convert COO adjacency list to CSR format"""
    # Create the CSR matrix
    col_index = []
    row_index = np.zeros(n_points + 1, dtype=int)

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

def rcm(problem: Problem) -> Problem:

    coo_adj = problem_coo(problem)
    csr_col, csr_row = csr_from_coo(coo_adj, problem.points.shape[0])

    # Temporarily replace with scipy implementation
    # from scipy.sparse import csr_matrix
    # from scipy.sparse.csgraph import reverse_cuthill_mckee

    # A = csr_matrix((np.ones(len(csr_col)), csr_col, csr_row))
    # order = reverse_cuthill_mckee(A)
    # order = list(order)



    min_deg = (None, np.inf)
    for i in range(len(csr_row) - 1):
        deg = csr_row[i + 1] - csr_row[i]
        if deg < min_deg[1]:
            min_deg = (i, deg)
    if min_deg[0] is None:
        # Return an error if no minimum degree node is found
        raise ValueError("No minimum degree node found in mesh, have fun debugging!")

    queue = []
    order = [min_deg[0]]
    visited = np.zeros(len(csr_row) - 1, dtype=bool)
    visited[min_deg[0]] = True
    for i in range(csr_row[min_deg[0]], csr_row[min_deg[0] + 1]):
        queue.append(csr_col[i])
    # Sort the queue by degree
    queue.sort(key=lambda x: csr_row[x + 1] - csr_row[x])
    while queue:
        node = queue.pop(0)
        if not visited[node]:
            order.append(node)
            visited[node] = True
            temp_queue = []
            for i in range(csr_row[node], csr_row[node + 1]):
                # if csr_col[i] not in order:
                if not visited[csr_col[i]]:
                    temp_queue.append(csr_col[i])
            temp_queue.sort(key=lambda x: csr_row[x + 1] - csr_row[x])
            queue = queue + temp_queue
    # Reverse the order to get the RCM ordering
    order.reverse()

    # inverse_order = np.zeros(len(order), dtype=int)
    # for i in range(len(order)):
    #     inverse_order[order[i]] = i

    # Renumber
    # for edge in problem.edges:
    #     for i in range(2):
    #         edge[i] = order.index(edge[i])
    # for edge in problem.edges:
    #     for i in range(2):
    #         edge[i] = inverse_order[edge[i]]

    # for element in problem.elements:
    #     for i in range(problem.n_local):
    # #         element[i] = order.index(element[i])
    # for element in problem.elements:
    #     for i in range(problem.n_local):
    #         element[i] = inverse_order[element[i]]

    # problem.points = problem.points[order]

    problem.order = order
    problem.renum = list(np.zeros(len(order), dtype=int))
    for i in range(len(order)):
        problem.renum[order[i]] = i

    return problem
