import numpy as np
import matplotlib.pyplot as plt

def plot_point_cloud_displacement(problem, solution, deformation_factor=1e5):
    """
    Plots the mesh as a point cloud of deformed node coordinates.
    Each node is colored according to its displacement magnitude.

    Parameters
    ----------
    problem : Problem
        An object that contains
            - points: (n_points, 2) array of original node coordinates
            - n_points: number of nodes
    solution : np.ndarray
        1D array of length 2*n_points containing [u_x0, u_y0, u_x1, u_y1, ...].
    deformation_factor : float
        Factor by which to scale the displacement for visualization.
    """
    # Number of points in the mesh
    n_points = problem.n_points

    # Reshape the solution so that U[i, 0] = x-displacement at node i
    # and U[i, 1] = y-displacement at node i
    U = solution.reshape((n_points, 2))

    # Copy original coordinates to avoid overwriting them
    deformed_points = problem.points.copy()

    # Scale the coordinates by the displacement
    # X[i] += u_x[i] * deformation_factor,  Y[i] += u_y[i] * deformation_factor
    deformed_points[:, 0] += U[:, 0] * deformation_factor
    deformed_points[:, 1] += U[:, 1] * deformation_factor

    # Displacement magnitude at each node
    displacement_magnitude = np.sqrt(U[:, 0]**2 + U[:, 1]**2)

    print("Maximum displacement magnitude: {:.10e}".format(max(displacement_magnitude)))

    # Plot
    plt.figure()  # separate figure for this plot
    sc = plt.scatter(
        deformed_points[:, 0],
        deformed_points[:, 1],
        c=displacement_magnitude,   # color by displacement magnitude
        cmap='jet',                 # choose any colormap you like
        marker='o'
    )
    plt.colorbar(sc, label='Displacement magnitude')
    plt.xlabel("x (deformed)")
    plt.ylabel("y (deformed)")
    plt.title("Deformed Mesh (scaled by {:.1e})".format(deformation_factor))
    plt.axis("equal")  # maintain aspect ratio
    plt.show()


def plot_point_cloud(problem) -> None:
    """
    Plots the mesh as a simple point cloud of all node coordinates.
    """
    plt.figure()  # each plot in its own figure
    plt.scatter(problem.points[:, 0], problem.points[:, 1], marker='o')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Simple point cloud of the mesh nodes")
    plt.axis("equal")  # so x and y scales are the same
    plt.show()


def plot_mesh(problem) -> None:
    """
    Plots the mesh by drawing each element.
    - For triangles, we connect nodes (n0->n1->n2->n0).
    - For quads, we connect nodes (n0->n1->n2->n3->n0).
    Also plots all the nodes as a scatter for clarity.
    """
    plt.figure()  # Make sure each chart is its own figure

    # First, scatter all node coordinates
    plt.scatter(problem.points[:, 0], problem.points[:, 1], marker='o')

    # Draw each element
    if problem.n_local == 3:
        # Triangles
        for elem_indices in problem.elements:
            n0, n1, n2 = elem_indices
            x0, y0 = problem.points[n0]
            x1, y1 = problem.points[n1]
            x2, y2 = problem.points[n2]
            # draw the lines n0->n1, n1->n2, n2->n0
            plt.plot([x0, x1], [y0, y1])
            plt.plot([x1, x2], [y1, y2])
            plt.plot([x2, x0], [y2, y0])
    elif problem.n_local == 4:
        # Quads
        for elem_indices in problem.elements:
            n0, n1, n2, n3 = elem_indices
            x0, y0 = problem.points[n0]
            x1, y1 = problem.points[n1]
            x2, y2 = problem.points[n2]
            x3, y3 = problem.points[n3]
            # draw the lines n0->n1, n1->n2, n2->n3, n3->n0
            plt.plot([x0, x1], [y0, y1])
            plt.plot([x1, x2], [y1, y2])
            plt.plot([x2, x3], [y2, y3])
            plt.plot([x3, x0], [y3, y0])

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Mesh visualization")
    plt.axis("equal")  # Keep aspect ratio so elements look correct
    plt.show()
