from problem import Problem
from parser import read_mesh

def main():
    problem = read_mesh("./data/elasticity.txt")

    E = 211e9
    nu = 0.3
    rho = 7.85e3
    g = 9.81

    iCase = "PLANAR_STRESS"
    if (iCase == "PLANAR_STRESS"):
        problem.A = E / (1 - nu * nu);
        problem.B = E * nu / (1 - nu * nu);
        problem.C = E / (2 * (1 + nu));
    elif (iCase == "PLANAR_STRAIN" or iCase == "AXISYM"):
        problem.A = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
        problem.B = E * nu / ((1 + nu) * (1 - 2 * nu));
        problem.C = E / (2 * (1 + nu));

    # addBoundaryConditions(problem, "Symmetry", "DIRICHLET_X", 0.0)
    # addBoundaryConditions(problem, "Symmetry", "DIRICHLET_Y", 0.0)

    # solution = solve_elasticity(problem)
