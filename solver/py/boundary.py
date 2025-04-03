from problem import Problem, Domain
import numpy as np

class BoundaryCondition:
    def __init__(self, value: float, type: str, domain: str):
        self.value = value
        self.type = type
        self.domain = domain

def add_boundary_condition(problem: Problem, domain: Domain, bc_type: str, value: float):
    problem.conditions.append(BoundaryCondition(value, bc_type, domain.name))

def apply_dirichlet_bc(A: np.ndarray, B: np.ndarray, dof: int, value: float) -> None:
    """
    Enforces U[dof] = value in a symmetric fashion:
      1) B[i] -= value*A[i,dof] for all i, then A[i,dof] = 0
      2) A[dof,j] = 0 for all j
      3) A[dof,dof] = 1
      4) B[dof] = value
    """
    size = A.shape[0]
    # (1) Update the entire B[] for the existing column
    # we know the value of U[dof] so we can already update B[i]
    # just as if U[dof] wasn't a variable
    for i in range(size):
        B[i] -= value * A[i, dof]
        A[i, dof] = 0.0

    # (2) Zero out row 'dof'
    for j in range(size):
        A[dof, j] = 0.0

    # (3) Set diagonal to 1
    A[dof, dof] = 1.0

    # (4) Set the new RHS
    B[dof] = value

def apply_dirichlet_bc_sym_banded(A: np.ndarray, B: np.ndarray, dof: int, value: float, bandwidth: int) -> None:
    """
    Enforces U[dof] = value in a symmetric fashion:
      1) B[i] -= value*A[i,dof] for all i, then A[i,dof] = 0
      2) A[dof,j] = 0 for all j
      3) A[dof,dof] = 1
      4) B[dof] = value
    """
    idx = lambda i, j, bw: (i + 1) * bw + j
    size = len(B)
    # (1) Update the entire B[] for the existing column
    # we know the value of U[dof] so we can already update B[i]
    # just as if U[dof] wasn't a variable
    # for i in range(size):
    #     B[i] -= value * A[i, dof]
    #     A[i, dof] = 0.0

    lower_bound = max(0, dof - bandwidth)
    upper_bound = min(size, dof + bandwidth + 1)
    for i in range(lower_bound, upper_bound):
        if dof > i:
            B[i] -= value * A[idx(dof, i, bandwidth)]
            A[idx(dof, i, bandwidth)] = 0.0
            continue
        B[i] -= value * A[idx(i, dof, bandwidth)]
        A[idx(i, dof, bandwidth)] = 0.0

    # (2) Zero out row 'dof'
    # for j in range(size):
    #     A[dof, j] = 0.0

    for j in range(lower_bound, dof):
        A[idx(dof, j, bandwidth)] = 0.0

    # (3) Set diagonal to 1
    # A[dof, dof] = 1.0

    A[idx(dof, dof, bandwidth)] = 1.0

    # (4) Set the new RHS
    B[dof] = value

def apply_dirichlet_bc_coo(coo: np.ndarray, B: np.ndarray, dof: int, value: float, nnz: int) -> None:
    """
    Enforces U[dof] = value in a symmetric fashion:
      1) B[i] -= value*A[i,dof] for all i, then A[i,dof] = 0
      2) A[dof,j] = 0 for all j
      3) A[dof,dof] = 1
      4) B[dof] = value
    """
    pass


def apply_neumann(problem, B):
    """
    Assembles the Neumann (force) contributions on boundary edges
    into the global RHS vector. This directly mirrors the C routine
    `femElasticityAssembleNeumann()`.
    """
    # Unpack needed attributes from `problem`
    B_global       = B                         # or however you store your system's RHS
    edge_rule      = problem.edge_integration_rule
    edge_shapes    = problem.edge_shapes       # shape functions for edges (2-node)
    edges          = problem.edges             # shape (n_edges, 2)
    points         = problem.points            # shape (n_points, 2)
    conditions     = problem.conditions        # list of boundary conditions
    n_local_edge   = 2  # for a 2-node edge

    # Loop over all boundary conditions
    for bc in conditions:
        bc_type  = bc.type
        bc_value = bc.value

        # Only handle NEUMANN cases; skip DIRICHLET
        if bc_type == "DIRICHLET_X" or bc_type == "DIRICHLET_Y":
            continue

        # Decide if it's NEUMANN_X or NEUMANN_Y
        if bc_type == "NEUMANN_X":
            shift = 0
        elif bc_type == "NEUMANN_Y":
            shift = 1
        else:
            # If you have other types like traction in normal/tangent,
            # you'd handle them differently. For now, skip unknown types.
            continue

        domain = [i for i in problem.domains if i.name == bc.domain][0]
        print(domain.name)

        for iEdge in domain.elem:
            # Extract the 2 nodes that define this edge
            # e.g. edges[iEdge] => [node0, node1]
            map_ = edges[iEdge]  # shape (2,)

            if problem.renum != []:
                map_ = [problem.renum[i] for i in map_]

            # Coordinates of the 2 edge nodes
            x_e = points[map_, 0]  # array of length 2
            y_e = points[map_, 1]  # array of length 2

            # Loop over the Gauss points on this edge
            for k in range(edge_rule.n_points):
                xsi    = edge_rule.points[k, 0]   # 1D gauss point
                weight = edge_rule.weights[k]

                # Evaluate shape functions & derivative wrt xsi
                phi   = edge_shapes.shape_fun(xsi)            # shape (2,)
                dphix = edge_shapes.shape_fun_derivatives(xsi)  # shape (2,)

                # Compute the Jacobian = length element in 1D
                # For a 2-node edge on [-1,1], dxdxsi = sum_i( x_e[i]* dphix[i] )
                dxdxsi = 0.0
                for i_loc in range(n_local_edge):
                    dxdxsi += x_e[i_loc] * dphix[i_loc]
                # The "length" = abs(dxdxsi)
                jac = abs(dxdxsi)

                # Accumulate into the global RHS
                for i_loc in range(n_local_edge):
                    i_node = map_[i_loc]
                    B_global[2*i_node + shift] += phi[i_loc] * bc_value * weight * jac

    return B_global
