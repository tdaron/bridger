from problem import Problem
from parser import read_mesh

from boundary import add_boundary_condition, apply_neumann
from solver import np_solve
from assemble import assemble_system

from plot import plot_point_cloud_displacement

def main():
    problem = read_mesh("./data/elasticity.txt")

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

    symmetry = [i for i in problem.domains if i.name == "Symmetry"][0]
    bottom = [i for i in problem.domains if i.name == "Bottom"][0]
    top = [i for i in problem.domains if i.name == "Top"][0]

    add_boundary_condition(problem, symmetry, "DIRICHLET_X", 0.0)
    add_boundary_condition(problem, bottom, "DIRICHLET_Y", 0.0)
    add_boundary_condition(problem, top, "NEUMANN_Y", -1e4);


    A, B = assemble_system(problem)

    import matplotlib.pyplot as plt
    plt.spy(A)
    plt.show()

    B = apply_neumann(problem, B)
    solution = np_solve(A, B)

    plot_point_cloud_displacement(problem, solution, deformation_factor=1e5)

if __name__ == "__main__":
    main()
