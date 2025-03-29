import numpy as np
from numpy.typing import NDArray
from typing import List, Tuple

class IntegrationRule:
    def __init__(self, n_points: int, points: NDArray[np.float64], weights: NDArray[np.float64]):
        """
        points: shape (n_points, dim)
          For quads/triangles, dim=2 (xsi, eta).
          For edges, dim=1 (xsi).
        weights: shape (n_points,)
        """
        self.n_points = n_points
        self.points = points
        self.weights = weights

    def __str__(self):
        return f"Integration rule with {self.n_points} points"

#
# -------------------- Integration Rules --------------------
#

def create_integration_rule_4point_quad() -> IntegrationRule:
    """
    2x2 Gauss rule on square [-1,1]x[-1,1].
    """
    # Each row is [xsi, eta]
    points = np.array([
        [-0.577350269189626,  0.577350269189626],
        [-0.577350269189626, -0.577350269189626],
        [ 0.577350269189626, -0.577350269189626],
        [ 0.577350269189626,  0.577350269189626],
    ], dtype=np.float64)
    weights = np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float64)
    return IntegrationRule(n_points=4, points=points, weights=weights)

def create_integration_rule_3point_triangle() -> IntegrationRule:
    """
    3-point rule on reference triangle with corners (0,0), (1,0), (0,1).
    """
    # Each row is [xsi, eta] in barycentric form
    points = np.array([
        [1.0 / 6.0, 1.0 / 6.0],
        [2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0],
    ], dtype=np.float64)
    weights = np.array([1.0 / 6.0]*3, dtype=np.float64)
    return IntegrationRule(n_points=3, points=points, weights=weights)

def create_integration_rule_2point_edge() -> IntegrationRule:
    """
    2-point Gauss-Legendre rule on [-1,1] for edges.
    """
    # Each row is just [xsi], so shape is (2,1)
    points = np.array([
        [ 0.577350269189626],
        [-0.577350269189626],
    ], dtype=np.float64)
    weights = np.array([1.0, 1.0], dtype=np.float64)
    return IntegrationRule(n_points=2, points=points, weights=weights)
