import numpy as np
from numpy.typing import NDArray
from typing import List

import matplotlib.pyplot as plt

from integration import create_integration_rule_3point_triangle, create_integration_rule_2point_edge, create_integration_rule_4point_quad
from shape import Shapes

class Domain:
    def __init__(self, name: str, elem: NDArray[np.int_]):
        self.name = name
        self.elem = elem
        self.n_elems = elem.shape[0]

class Problem:
    def __init__(
        self,
        points: NDArray[np.float64],
        n_points: int,
        edges: NDArray[np.int_],
        n_edges: int,
        elements: NDArray[np.int_],
        n_elements: int,
        n_local: int,
        domains: List[Domain]
    ):
        self.points = points         # shape: (n_points, 2)
        self.n_points = n_points
        self.edges = edges           # shape: (n_edges, 2)
        self.n_edges = n_edges
        self.elements = elements     # shape: (n_elements, n_local)
        self.n_elements = n_elements
        self.n_local = n_local
        self.domains = domains

        self.element_shapes = Shapes(n_local)
        self.edge_shapes = Shapes(2)

        if n_local == 3:
            self.elem_integration_rule = create_integration_rule_3point_triangle()
        elif n_local == 4:
            self.elem_integration_rule = create_integration_rule_4point_quad()
        else:
            raise ValueError(f"Unknown element type: {n_local}")

        self.edge_integration_rule = create_integration_rule_2point_edge()

        self.A, self.B, self.C = 0.0, 0.0, 0.0

        self.node_constraints = np.ones(n_points * 2, dtype=np.int_) * -1
        self.conditions = []

        self.bandwidth = np.inf
        self.renum = []
        self.order = []

    def __str__(self):
        name = "triangles" if self.n_local == 3 else "quads"
        return (
            f"Problem with {self.n_points} points, {self.n_edges} edges, "
            f"{self.n_elements} {name}, and {len(self.domains)} domains"
        )

    def print_geometry(self) -> None:
        # 1) Print Node (Point) information
        print(f"Number of nodes {self.n_points}")
        for i in range(self.n_points):
            x, y = self.points[i]
            print(f"{i:6d} : {x:14.7e} {y:14.7e}")

        # 2) Print Edge information
        print(f"Number of edges {self.n_edges}")
        for i in range(self.n_edges):
            n0, n1 = self.edges[i]
            print(f"{i:6d} : {n0:6d} {n1:6d}")

        # 3) Print Element information (triangles vs quads)
        if self.n_local == 3:
            print(f"Number of triangles {self.n_elements}")
            for i in range(self.n_elements):
                n0, n1, n2 = self.elements[i]
                print(f"{i:6d} : {n0:6d} {n1:6d} {n2:6d}")
        elif self.n_local == 4:
            print(f"Number of quads {self.n_elements}")
            for i in range(self.n_elements):
                n0, n1, n2, n3 = self.elements[i]
                print(f"{i:6d} : {n0:6d} {n1:6d} {n2:6d} {n3:6d}")

        # 4) Print Domain information
        print(f"Number of domains {len(self.domains)}")
        for i_domain, domain in enumerate(self.domains):
            print(f"  Domain : {i_domain:6d}")
            print(f"  Name : {domain.name}")
            print(f"  Number of elements : {domain.n_elems:6d}")

            # Print the domain's element indices, 10 per line
            for i, elem_idx in enumerate(domain.elem):
                print(f"{elem_idx:6d}", end="")
                if (i + 1) != domain.n_elems and (i + 1) % 10 == 0:
                    print()  # line break after 10 entries
            print("\n", end="")  # final line break after each domain
