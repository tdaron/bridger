from problem import Problem

import numpy as np
from gmsh.gmsh-4.13.1-MacOSARM-sdk.share.doc.gmsh.tutorials.python.x2 import coords

# x axis renumbering
def renumber_x_axis(problem: Problem) -> Problem:
    id_coords_list = list(enumerate(problem.points))

    # sort by x coordinate
    id_coords_list.sort(key=lambda x: x[1][0])

    # renumber
    for edge in problem.edges:
        for i in range(2):
            edge[i] = [x[0] for x in id_coords_list].index(edge[i])

    for element in problem.elements:
        for i in range(problem.n_local):
            element[i] = [x[0] for x in id_coords_list].index(element[i])

    problem.points = np.array([x[1] for x in id_coords_list])

    return problem

def rcm(problem: Problem) -> Problem:
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

    prev_row = coo_adj[0][0]
    prev = None
    for i in range(len(coo_adj)):
        # Skip duplicates
        if prev == coo_adj[i]:
            continue
        col_index.append(coo_adj[i][1])
        # Set all absent row indices to the current column index
        for j in range(prev_row + 1, coo_adj[i][0]):
            row_index[j] = len(col_index) - 1
        prev_row = coo_adj[i][0]
        prev = coo_adj[i]
    row_index[-1] = len(col_index)
