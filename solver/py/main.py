from problem import Problem
from parse import read_mesh

from boundary import add_boundary_condition, apply_neumann
from solver import np_solve
from assemble import assemble_system

from band import renumber_x_axis

from plot import plot_point_cloud_displacement

def main():
    problem = read_mesh("./data/mesh.txt")

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

    # symmetry = [i for i in problem.domains if i.name == "Symmetry"][0]
    # bottom = [i for i in problem.domains if i.name == "Bottom"][0]
    # top = [i for i in problem.domains if i.name == "Top"][0]

    # add_boundary_condition(problem, symmetry, "DIRICHLET_X", 0.0)
    # add_boundary_condition(problem, bottom, "DIRICHLET_Y", 0.0)
    # add_boundary_condition(problem, top, "NEUMANN_Y", -1e4);

    pillar1 = [i for i in problem.domains if i.name == "Pillar1"][0]
    pillar2 = [i for i in problem.domains if i.name == "Pillar2"][0]
    pillar3 = [i for i in problem.domains if i.name == "Pillar3"][0]
    pillar4 = [i for i in problem.domains if i.name == "Pillar4"][0]

    left = [i for i in problem.domains if i.name == "Left"][0]
    right = [i for i in problem.domains if i.name == "Right"][0]

    add_boundary_condition(problem, pillar1, "DIRICHLET_Y", 0.0)
    add_boundary_condition(problem, pillar2, "DIRICHLET_Y", 0.0)
    add_boundary_condition(problem, pillar3, "DIRICHLET_Y", 0.0)
    add_boundary_condition(problem, pillar4, "DIRICHLET_Y", 0.0)

    add_boundary_condition(problem, left, "DIRICHLET_X", 0.0)
    add_boundary_condition(problem, right, "DIRICHLET_X", 0.0)


    A, B = assemble_system(problem)

    # calculate band width
    bandwidth = 0
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if A[i, j] != 0:
                bandwidth = max(bandwidth, abs(i - j))
    print(f"Bandwidth: {bandwidth}")

    import matplotlib.pyplot as plt
    plt.spy(A)
    plt.show()

    # x axis renumbering
    problem = renumber_x_axis(problem)

    A, B = assemble_system(problem)

    # calculate band width
    bandwidth = 0
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if A[i, j] != 0:
                bandwidth = max(bandwidth, abs(i - j))
    print(f"Bandwidth: {bandwidth}")

    plt.spy(A)
    plt.show()

    B = apply_neumann(problem, B)
    solution = np_solve(A, B)

    plot_point_cloud_displacement(problem, solution, deformation_factor=5e3)

if __name__ == "__main__":
    main()
