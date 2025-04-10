import numpy as np
from numpy.typing import NDArray

from typing import Tuple, List

from boundary import apply_dirichlet_bc, apply_dirichlet_bc_sym_banded

god_damn = np.zeros((1, 1), dtype=np.float64)

def assemble_system(problem) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    n_points = problem.n_points
    n_elems  = problem.n_elements
    n_local  = problem.n_local
    n_dof    = 2 * n_points

    # Get Material Parameters
    a = problem.A
    b = problem.B
    c = problem.C

    # Gravity Parameters
    rho = 7.85e3
    g = 9.81

    A_global = np.zeros((n_dof, n_dof), dtype=np.float64)
    B_global = np.zeros(n_dof, dtype=np.float64)

    # Shorthands for shape function calls
    shape_fun            = problem.element_shapes.shape_fun
    shape_fun_deriv      = problem.element_shapes.shape_fun_derivatives
    integration_rule     = problem.elem_integration_rule

    for e in range(n_elems):
        # 'map' will hold the global node indices for this element
        map_ = problem.elements[e]

        x_e = problem.points[map_, 0]
        y_e = problem.points[map_, 1]

        if problem.renum != []:
            map_ = [problem.renum[i] for i in map_]

        # Loop over the integration points in the element
        for k in range(integration_rule.n_points):
            xsi    = integration_rule.points[k, 0]
            eta    = integration_rule.points[k, 1]
            weight = integration_rule.weights[k]

            # Evaluate shape functions & their ref-derivatives at (xsi, eta)
            phi = shape_fun(xsi, eta)  # length n_local
            dphidxsi, dphideta = shape_fun_deriv(xsi, eta)

            # Compute the Jacobian (xsi->X, eta->Y) to get dxdxsi, dxdeta, ...
            dxdxsi = 0.0; dxdeta = 0.0
            dydxsi = 0.0; dydeta = 0.0

            for i_loc in range(n_local):
                dxdxsi += x_e[i_loc] * dphidxsi[i_loc]
                dxdeta += x_e[i_loc] * dphideta[i_loc]
                dydxsi += y_e[i_loc] * dphidxsi[i_loc]
                dydeta += y_e[i_loc] * dphideta[i_loc]

            J_e = abs(dxdxsi * dydeta - dxdeta * dydxsi)

            dphidx = np.zeros(n_local)
            dphidy = np.zeros(n_local)

            for i_loc in range(n_local):
                dphidx[i_loc] = ( dphidxsi[i_loc] * dydeta - dphideta[i_loc] * dydxsi ) / J_e
                dphidy[i_loc] = ( dphideta[i_loc] * dxdxsi - dphidxsi[i_loc] * dxdeta ) / J_e

            for i_loc in range(n_local):
                i_node = map_[i_loc]
                i_dof_x = 2 * i_node
                i_dof_y = 2 * i_node + 1

                for j_loc in range(n_local):
                    j_node = map_[j_loc]
                    j_dof_x = 2 * j_node
                    j_dof_y = 2 * j_node + 1

                    A_global[i_dof_x, j_dof_x] += (a * dphidx[i_loc] * dphidx[j_loc]
                                                   + c * dphidy[i_loc] * dphidy[j_loc]) * weight * J_e
                    A_global[i_dof_x, j_dof_y] += (b * dphidx[i_loc] * dphidy[j_loc]
                                                   + c * dphidy[i_loc] * dphidx[j_loc]) * weight * J_e
                    A_global[i_dof_y, j_dof_x] += (c * dphidx[i_loc] * dphidy[j_loc]
                                                   + b * dphidy[i_loc] * dphidx[j_loc]) * weight * J_e
                    A_global[i_dof_y, j_dof_y] += (c * dphidx[i_loc] * dphidx[j_loc]
                                                   + a * dphidy[i_loc] * dphidy[j_loc]) * weight * J_e
                B_global[i_dof_x] += 0.0
                B_global[i_dof_y] -= (rho * g * phi[i_loc] * weight * J_e)

    # Count and print the number of non-zero elements in the global matrix within a tolerance
    print("Number of non-zero elements in the global matrix: ", np.count_nonzero(A_global > 1e-6))

    global god_damn
    god_damn = A_global.copy()

    # Apply Dirichlet BCs
    for bc in problem.conditions:
        if bc.type not in ["DIRICHLET_X", "DIRICHLET_Y"]:
            continue
        domain = [i for i in problem.domains if i.name == bc.domain][0]
        for iEdge in domain.elem:
            map_ = problem.edges[iEdge]
            if problem.renum != []:
                map_ = [problem.renum[i] for i in map_]
            shift = 0 if bc.type == "DIRICHLET_X" else 1
            for i_loc in range(2):
                i_node = map_[i_loc]
                i_dof_x = 2 * i_node
                i_dof_y = 2 * i_node + 1
                apply_dirichlet_bc(A_global, B_global, i_dof_x + shift, bc.value)


    return A_global, B_global


def assemble_system_banded(problem) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    idx = lambda i, j: (i + 1) * problem.bandwidth + j


    n_points = problem.n_points
    n_elems  = problem.n_elements
    n_local  = problem.n_local
    n_dof    = 2 * n_points

    # Get Material Parameters
    a = problem.A
    b = problem.B
    c = problem.C

    # Gravity Parameters
    rho = 7.85e3
    g = 9.81

    # A_global = np.zeros((n_dof, n_dof), dtype=np.float64)
    # A_global = np.zeros((n_dof, problem.bandwidth * 2 + 1), dtype=np.float64)
    A_global = np.zeros(n_dof * (problem.bandwidth + 1), dtype=np.float64)
    B_global = np.zeros(n_dof, dtype=np.float64)

    # Shorthands for shape function calls
    shape_fun            = problem.element_shapes.shape_fun
    shape_fun_deriv      = problem.element_shapes.shape_fun_derivatives
    integration_rule     = problem.elem_integration_rule

    for e in range(n_elems):
        # 'map' will hold the global node indices for this element
        map_ = problem.elements[e]

        x_e = problem.points[map_, 0]
        y_e = problem.points[map_, 1]

        if problem.renum != []:
                    map_ = [problem.renum[i] for i in map_]


        # Loop over the integration points in the element
        for k in range(integration_rule.n_points):
            xsi    = integration_rule.points[k, 0]
            eta    = integration_rule.points[k, 1]
            weight = integration_rule.weights[k]

            # Evaluate shape functions & their ref-derivatives at (xsi, eta)
            phi = shape_fun(xsi, eta)  # length n_local
            dphidxsi, dphideta = shape_fun_deriv(xsi, eta)

            # Compute the Jacobian (xsi->X, eta->Y) to get dxdxsi, dxdeta, ...
            dxdxsi = 0.0; dxdeta = 0.0
            dydxsi = 0.0; dydeta = 0.0

            for i_loc in range(n_local):
                dxdxsi += x_e[i_loc] * dphidxsi[i_loc]
                dxdeta += x_e[i_loc] * dphideta[i_loc]
                dydxsi += y_e[i_loc] * dphidxsi[i_loc]
                dydeta += y_e[i_loc] * dphideta[i_loc]

            J_e = abs(dxdxsi * dydeta - dxdeta * dydxsi)

            dphidx = np.zeros(n_local)
            dphidy = np.zeros(n_local)

            for i_loc in range(n_local):
                dphidx[i_loc] = ( dphidxsi[i_loc] * dydeta - dphideta[i_loc] * dydxsi ) / J_e
                dphidy[i_loc] = ( dphideta[i_loc] * dxdxsi - dphidxsi[i_loc] * dxdeta ) / J_e

            for i_loc in range(n_local):
                i_node = map_[i_loc]
                i_dof_x = 2 * i_node
                i_dof_y = 2 * i_node + 1

                for j_loc in range(n_local):
                    j_node = map_[j_loc]
                    j_dof_x = 2 * j_node
                    j_dof_y = 2 * j_node + 1

                    if j_dof_x <= i_dof_x:
                        A_global[idx(i_dof_x, j_dof_x)] += (a * dphidx[i_loc] * dphidx[j_loc]
                                                   + c * dphidy[i_loc] * dphidy[j_loc]) * weight * J_e
                    if j_dof_y <= i_dof_x:
                        A_global[idx(i_dof_x, j_dof_y)] += (b * dphidx[i_loc] * dphidy[j_loc]
                                                   + c * dphidy[i_loc] * dphidx[j_loc]) * weight * J_e
                    if j_dof_x <= i_dof_y:
                        A_global[idx(i_dof_y, j_dof_x)] += (c * dphidx[i_loc] * dphidy[j_loc]
                                                   + b * dphidy[i_loc] * dphidx[j_loc]) * weight * J_e
                    if j_dof_y <= i_dof_y:
                        A_global[idx(i_dof_y, j_dof_y)] += (c * dphidx[i_loc] * dphidx[j_loc]
                                                    + a * dphidy[i_loc] * dphidy[j_loc]) * weight * J_e
                    B_global[i_dof_x] += 0.0
                B_global[i_dof_y] -= (rho * g * phi[i_loc] * weight * J_e)

    # Apply Dirichlet BCs
    for bc in problem.conditions:
        if bc.type not in ["DIRICHLET_X", "DIRICHLET_Y"]:
            continue
        domain = [i for i in problem.domains if i.name == bc.domain][0]
        for iEdge in domain.elem:
            map_ = problem.edges[iEdge]
            if problem.renum != []:
                map_ = [problem.renum[i] for i in map_]
            shift = 0 if bc.type == "DIRICHLET_X" else 1
            for i_loc in range(2):
                i_node = map_[i_loc]
                i_dof_x = 2 * i_node
                i_dof_y = 2 * i_node + 1
                apply_dirichlet_bc_sym_banded(A_global, B_global, i_dof_x + shift, bc.value, problem.bandwidth)


    return A_global, B_global

def assemble_csr_system(problem):
    n_points = problem.n_points
    n_elems  = problem.n_elements
    n_local  = problem.n_local
    n_dof    = 2 * n_points

    # Get Material Parameters
    a = problem.A
    b = problem.B
    c = problem.C

    # Gravity Parameters
    rho = 7.85e3
    g = 9.81

    # A_global = np.zeros((n_dof, n_dof), dtype=np.float64)
    coo_global = np.zeros((n_elems * n_local * n_local * 4, 3), dtype=np.int64)
    coo_vals = np.zeros((n_elems * n_local * n_local * 4), dtype=np.float64)
    B_global = np.zeros(n_dof, dtype=np.float64)

    # Shorthands for shape function calls
    shape_fun            = problem.element_shapes.shape_fun
    shape_fun_deriv      = problem.element_shapes.shape_fun_derivatives
    integration_rule     = problem.elem_integration_rule

    nnz = 0
    for e in range(n_elems):
        # 'map' will hold the global node indices for this element
        map_ = problem.elements[e]

        x_e = problem.points[map_, 0]
        y_e = problem.points[map_, 1]

        if problem.renum != []:
            map_ = [problem.renum[i] for i in map_]

        local_nnz = 0
        # Loop over the integration points in the element
        for k in range(integration_rule.n_points):
            xsi    = integration_rule.points[k, 0]
            eta    = integration_rule.points[k, 1]
            weight = integration_rule.weights[k]

            # Evaluate shape functions & their ref-derivatives at (xsi, eta)
            phi = shape_fun(xsi, eta)  # length n_local
            dphidxsi, dphideta = shape_fun_deriv(xsi, eta)

            # Compute the Jacobian (xsi->X, eta->Y) to get dxdxsi, dxdeta, ...
            dxdxsi = 0.0; dxdeta = 0.0
            dydxsi = 0.0; dydeta = 0.0

            for i_loc in range(n_local):
                dxdxsi += x_e[i_loc] * dphidxsi[i_loc]
                dxdeta += x_e[i_loc] * dphideta[i_loc]
                dydxsi += y_e[i_loc] * dphidxsi[i_loc]
                dydeta += y_e[i_loc] * dphideta[i_loc]

            J_e = abs(dxdxsi * dydeta - dxdeta * dydxsi)

            dphidx = np.zeros(n_local)
            dphidy = np.zeros(n_local)

            local_nnz = 0
            for i_loc in range(n_local):
                dphidx[i_loc] = ( dphidxsi[i_loc] * dydeta - dphideta[i_loc] * dydxsi ) / J_e
                dphidy[i_loc] = ( dphideta[i_loc] * dxdxsi - dphidxsi[i_loc] * dxdeta ) / J_e

            for i_loc in range(n_local):
                i_node = map_[i_loc]
                i_dof_x = 2 * i_node
                i_dof_y = 2 * i_node + 1

                for j_loc in range(n_local):
                    j_node = map_[j_loc]
                    j_dof_x = 2 * j_node
                    j_dof_y = 2 * j_node + 1

                    val_xx = (a * dphidx[i_loc] * dphidx[j_loc] + c * dphidy[i_loc] * dphidy[j_loc]) * weight * J_e
                    val_xy = (b * dphidx[i_loc] * dphidy[j_loc] + c * dphidy[i_loc] * dphidx[j_loc]) * weight * J_e
                    val_yx = (c * dphidx[i_loc] * dphidy[j_loc] + b * dphidy[i_loc] * dphidx[j_loc]) * weight * J_e
                    val_yy = (c * dphidx[i_loc] * dphidx[j_loc] + a * dphidy[i_loc] * dphidy[j_loc]) * weight * J_e

                    # if abs(val_xx) > 1e-12:
                    coo_global[nnz + local_nnz, 0] = i_dof_x
                    coo_global[nnz + local_nnz, 1] = j_dof_x
                    coo_global[nnz + local_nnz, 2] = val_xx
                    coo_vals[nnz + local_nnz] += val_xx
                    local_nnz += 1

                    # if abs(val_xy) > 1e-12:
                    coo_global[nnz + local_nnz, 0] = i_dof_x
                    coo_global[nnz + local_nnz, 1] = j_dof_y
                    coo_global[nnz + local_nnz, 2] = val_xy
                    coo_vals[nnz + local_nnz] += val_xy
                    local_nnz += 1

                    # if abs(val_yx) > 1e-12:
                    coo_global[nnz + local_nnz, 0] = i_dof_y
                    coo_global[nnz + local_nnz, 1] = j_dof_x
                    coo_global[nnz + local_nnz, 2] = val_yx
                    coo_vals[nnz + local_nnz] += val_yx
                    local_nnz += 1

                    # if abs(val_yy) > 1e-12:
                    coo_global[nnz + local_nnz, 0] = i_dof_y
                    coo_global[nnz + local_nnz, 1] = j_dof_y
                    coo_global[nnz + local_nnz, 2] = val_yy
                    coo_vals[nnz + local_nnz] += val_yy
                    local_nnz += 1
                B_global[i_dof_x] += 0.0
                B_global[i_dof_y] -= (rho * g * phi[i_loc] * weight * J_e)
        nnz += local_nnz
    # print("COO nnz = ", nnz)
    coo_global = coo_global[:nnz, :]
    coo_vals = coo_vals[:nnz]

    # To sort by first column, then second column
    sort_indices = np.lexsort((coo_global[:,1], coo_global[:,0]))

    # Apply the same sort indices to both arrays
    coo_global = coo_global[sort_indices]
    coo_vals = coo_vals[sort_indices]

    col_index, row_index, values = csr_from_coo(coo_global, coo_vals, nnz, n_points)
    print(np.count_nonzero(values))
    print(np.count_nonzero(abs(values) > 1e-6))

    A = np.zeros((n_dof, n_dof), dtype=np.float64)
    for i in range(len(row_index) - 1):
        for j in range(row_index[i], row_index[i + 1]):
            A[i, col_index[j]] = values[j]
    coo_global = A

    # Apply Dirichlet BCs
    for bc in problem.conditions:
        if bc.type not in ["DIRICHLET_X", "DIRICHLET_Y"]:
            continue
        domain = [i for i in problem.domains if i.name == bc.domain][0]
        for iEdge in domain.elem:
            map_ = problem.edges[iEdge]
            if problem.renum != []:
                map_ = [problem.renum[i] for i in map_]
            shift = 0 if bc.type == "DIRICHLET_X" else 1
            for i_loc in range(2):
                i_node = map_[i_loc]
                i_dof_x = 2 * i_node
                i_dof_y = 2 * i_node + 1
                apply_dirichlet_bc(coo_global, B_global, i_dof_x + shift, bc.value)


    return coo_global, B_global



def csr_from_coo(coo_adj: NDArray, coo_vals, nnz: int, n_points: int) -> tuple:
    # Create the CSR matrix
    col_index = np.zeros(nnz, dtype=int)
    row_index = np.zeros(n_points * 2 + 1, dtype=int)
    values = np.zeros(nnz, dtype=float)

    prev_row = 0
    prev = np.array([-1, -1], dtype=int)
    nnz_new = 0
    for i in range(len(coo_adj)):
        # Add values of duplicate entries
        # to the previous entry which is the first
        # of the duplicates and was added to the lists.
        # The is because we assume the input coo
        # lists are sorted by column then row.
        if np.array_equal(prev, coo_adj[i][:2]):
            # values[nnz_new - 1] += coo_adj[i][2]
            values[nnz_new - 1] += coo_vals[i]
            continue
        if abs(coo_vals[i]) > 1e-10:
            col_index[nnz_new] = coo_adj[i][1]
            values[nnz_new] = coo_vals[i]
            nnz_new += 1
            # Set all absent rows that were skipped to the
            # current column index, this way they don't
            # have any values associated with them.
            for j in range(prev_row, coo_adj[i][0]):
                row_index[j + 1] = nnz_new - 1
            prev_row = coo_adj[i][0]
            prev = coo_adj[i][:2].copy()


    # Make sure to update the last rows that might be empty
    for j in range(prev_row + 1, len(row_index)):
        row_index[j] = nnz_new

    # Remove unused elements
    col_index = col_index[:nnz_new]
    values = values[:nnz_new]

    return col_index, row_index, values
