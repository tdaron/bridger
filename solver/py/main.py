from problem import Problem
from parse import read_mesh

from boundary import add_boundary_condition, apply_neumann
from solver import np_solve, sym_band_cholesky, c_grad_full, pre_c_grad_full, incomplete_cholesky
from assemble import assemble_system, assemble_system_banded, assemble_csr_system

from band import renumber_x_axis, rcm
import scipy.sparse.linalg

from plot import plot_point_cloud_displacement

import numpy as np

import time

def extract_full_from_sym_banded(A, n, band):
    A_full = np.zeros((n, n))
    idx = lambda i, j: i * (band + 1) + j - i + band
    for i in range(n):
        for j in range(n):
            if abs(i - j) > band:
                continue
            if j > i:
                A_full[i, j] = A[idx(j, i)]
            else:
                A_full[i, j] = A[idx(i, j)]
    return A_full

def main():
    # problem = read_mesh("./data/elasticity.txt")
    problem = read_mesh("./data/mensh.txt")

    E = 211.e9
    nu = 0.3

    iCase = "PLANAR_STRAIN"
    if (iCase == "PLANAR_STRESS"):
        problem.A = E / (1.0 - nu * nu);
        problem.B = E * nu / (1.0 - nu * nu);
        problem.C = E / (2.0 * (1.0 + nu));
    elif (iCase == "PLANAR_STRAIN" or iCase == "AXISYM"):
        problem.A = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        problem.B = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        problem.C = E / (2.0 * (1.0 + nu));

    # pillar1 = [i for i in problem.domains if i.name == "Pillar1"][0]
    # pillar2 = [i for i in problem.domains if i.name == "Pillar2"][0]
    # pillar3 = [i for i in problem.domains if i.name == "Pillar3"][0]
    # pillar4 = [i for i in problem.domains if i.name == "Pillar4"][0]

    # left = [i for i in problem.domains if i.name == "Left"][0]
    # right = [i for i in problem.domains if i.name == "Right"][0]

    # add_boundary_condition(problem, pillar1, "DIRICHLET_Y", 0.0)
    # add_boundary_condition(problem, pillar2, "DIRICHLET_Y", 0.0)
    # add_boundary_condition(problem, pillar3, "DIRICHLET_Y", 0.0)
    # add_boundary_condition(problem, pillar4, "DIRICHLET_Y", 0.0)

    # add_boundary_condition(problem, left, "DIRICHLET_X", 0.0)
    # add_boundary_condition(problem, right, "DIRICHLET_X", 0.0)


    # symmetry = [i for i in problem.domains if i.name == "Symmetry"][0]
    # bottom = [i for i in problem.domains if i.name == "Bottom"][0]

    # add_boundary_condition(problem, symmetry, "DIRICHLET_X", 0.0)
    # add_boundary_condition(problem, bottom, "DIRICHLET_Y", 0.0)


    pillar1 = [i for i in problem.domains if i.name == "Pillar1"][0]
    pillar2 = [i for i in problem.domains if i.name == "Pillar2"][0]
    pillar3 = [i for i in problem.domains if i.name == "Pillar3"][0]
    pillar4 = [i for i in problem.domains if i.name == "Pillar4"][0]

    left = [i for i in problem.domains if i.name == "Extremity0"][0]
    right = [i for i in problem.domains if i.name == "Extremity1"][0]

    add_boundary_condition(problem, pillar1, "DIRICHLET_Y", 0.0)
    add_boundary_condition(problem, pillar2, "DIRICHLET_Y", 0.0)
    add_boundary_condition(problem, pillar3, "DIRICHLET_Y", 0.0)
    add_boundary_condition(problem, pillar4, "DIRICHLET_Y", 0.0)

    add_boundary_condition(problem, pillar1, "DIRICHLET_X", 0.0)
    add_boundary_condition(problem, pillar2, "DIRICHLET_X", 0.0)
    add_boundary_condition(problem, pillar3, "DIRICHLET_X", 0.0)
    add_boundary_condition(problem, pillar4, "DIRICHLET_X", 0.0)

    add_boundary_condition(problem, left, "DIRICHLET_X", 0.0)
    add_boundary_condition(problem, right, "DIRICHLET_X", 0.0)

    add_boundary_condition(problem, left, "DIRICHLET_Y", 0.0)
    add_boundary_condition(problem, right, "DIRICHLET_Y", 0.0)



    # RCM renumbering
    problem = rcm(problem)

    # Calculate bandwidth
    dist = 0
    for i in range(problem.n_elements):
        for j in range(problem.n_local):
            for k in range(problem.n_local):
                diff = abs(problem.renum[problem.elements[i][j]] - problem.renum[problem.elements[i][k]])
                dist = max(dist, diff)
    # Need to double the element id distance because of the 2 dofs per node
    # and add 1 for the diagonal which goes from a single element when there
    # is a single degree of freedom per node to a 2x2 block when there are 2.
    # For example, the matrix will always be at least of bandwidth 1.
    problem.bandwidth = 2 * dist + 1
    print(problem.bandwidth)

    A, B = assemble_system(problem)
    # exit(1)
    A, B = assemble_csr_system(problem)

    # # Print minimum eigenvalue
    # eigvals = np.linalg.eigvals(A)
    # print(f"Min eigenvalue: {min(eigvals)}")

    # A, B = assemble_system_banded(problem)
    # A_full = extract_full_from_sym_banded(A, len(B), bandwidth)
    # Make sure the eigenvalues are positive
    # eigvals = np.linalg.eigvals(A_full)
    # print(f"Min eigenvalue: {min(eigvals)}")

    B = apply_neumann(problem, B)

    # import matplotlib.pyplot as plt
    # plt.spy(A)
    # plt.show()

    t0 = time.time()

    # solution = sym_band_cholesky(A, B, np.zeros_like(B), problem.bandwidth)
    # solution = np_solve(A, B)
    # Calculte condition number
    # print(f"Condition number: {np.linalg.cond(A)}")
    # Condition number with jacobi preconditioning
    # print(f"Condition number: {np.linalg.cond(A @ np.diag(1 / np.diag(A)))}")
    # Calculte incomplete cholesky preconditioner
    # L = incomplete_cholesky(A)
    # Make non-zero entries of L that are zero in A zero
    # Calculte condition number
    # print(f"Condition number: {np.linalg.cond(np.linalg.inv(L @ L.T) @ A)}")
    # exit(1)

    solution = np_solve(A, B)
    # solution = pre_c_grad_full(A, B, np.ones_like(B), 1e-12)
    # solution = c_grad_full(A, B, np.zeros_like(B), 1e-6)
    # ret = scipy.sparse.linalg.cg(A, B, x0=np.zeros_like(B), rtol=1e-6)
    # print(ret[1])
    # solution = scipy.sparse.linalg.cg(A, B, x0=np.zeros_like(B), rtol=1e-6)[0]

    t1 = time.time()
    print(f"Time: {t1 - t0}")

    solution_reordered = np.zeros_like(solution)
    for i in range(len(solution)):
        solution_reordered[2 * problem.order[i // 2] + i % 2] = solution[i]
    solution = solution_reordered


    plot_point_cloud_displacement(problem, solution, deformation_factor=1e5)
    # plot_point_cloud_displacement(problem, solution, deformation_factor=5e3)

if __name__ == "__main__":
    main()
