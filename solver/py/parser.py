from problem import Domain, Problem

import numpy as np
from typing import List

from plotter import plot_point_cloud

# Equivalent function to the C version, using numpy arrays.
def read_mesh(filename: str) -> Problem:
    with open(filename, 'r') as f:
        # ---- Read Nodes ----
        # Expecting a line: "Number of nodes <n>"
        line = f.readline().strip()
        n_nodes = int(line.split()[-1])

        X = np.empty(n_nodes, dtype=np.float64)
        Y = np.empty(n_nodes, dtype=np.float64)
        for i in range(n_nodes):
            # Expecting lines like: "<node_id> : <x> <y>"
            line = f.readline().strip()
            parts = line.split()
            # parts[0] is node id, parts[1] is ":", parts[2] is x, parts[3] is y
            X[i] = float(parts[2])
            Y[i] = float(parts[3])
        points = np.column_stack((X, Y))

        # ---- Read Edges ----
        # Expecting a line: "Number of edges <n>"
        line = f.readline().strip()
        n_edges = int(line.split()[-1])
        edges = np.empty((n_edges, 2), dtype=np.int_)
        for i in range(n_edges):
            # Expecting lines like: "<edge_id> : <node1> <node2>"
            line = f.readline().strip()
            parts = line.split()
            edges[i, 0] = int(parts[2])
            edges[i, 1] = int(parts[3])

        # ---- Read Elements ----
        # Expecting a line: "Number of <element_type> <n>"
        line = f.readline().strip()
        parts = line.split()
        # For a line like: "Number of triangles 10"
        element_type = parts[2].lower()  # e.g. "triangles" or "quads"
        n_elements = int(parts[3])
        if element_type.startswith("tri"):
            n_local = 3
        elif element_type.startswith("quad"):
            n_local = 4
        else:
            raise ValueError(f"Unknown element type: {element_type}")
        elements = np.empty((n_elements, n_local), dtype=np.int_)
        for i in range(n_elements):
            # Expecting lines like: "<elem_id> : <n1> <n2> <n3>" (or 4 nodes for quads)
            line = f.readline().strip()
            parts = line.split()
            for j in range(n_local):
                elements[i, j] = int(parts[j+2])

        # ---- Read Domains ----
        # Expecting a line: "Number of domains <n>"
        line = f.readline().strip()
        n_domains = int(line.split()[-1])
        domains: List[Domain] = []
        for _ in range(n_domains):
            # Read domain header: "Domain : <id>" (we ignore the id)
            f.readline()
            # Next line: "Name : <name>"
            name_line = f.readline().strip()
            name = name_line.split(":", 1)[1].strip()
            # Next line: "Number of elements : <n>"
            n_line = f.readline().strip()
            n_dom_elems = int(n_line.split()[-1])

            # Read domain element indices.
            # They may span multiple lines. We continue reading until we have n_dom_elems numbers.
            elems = []
            while len(elems) < n_dom_elems:
                line = f.readline().strip()
                if not line:
                    continue
                elems.extend(int(x) for x in line.split())
            # Convert to numpy array
            domain_elem = np.array(elems[:n_dom_elems], dtype=np.int_)
            domains.append(Domain(name, domain_elem))

    return Problem(
        points,
        n_nodes,
        edges,
        n_edges,
        elements,
        n_elements,
        n_local,
        domains
    )

if __name__ == "__main__":
    problem = read_mesh("./data/elasticity.txt")
    print(problem)
    problem.print_geometry()
    problem.plot_point_cloud()
