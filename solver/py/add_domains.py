import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import sys
from typing import List, Dict, Tuple

class Domain:
    def __init__(self, name: str, elements: np.ndarray):
        self.name = name
        self.elements = elements

    def __repr__(self):
        return f"Domain(name='{self.name}', elements={self.elements})"

class Mesh:
    def __init__(self):
        self.nodes = None  # (n_nodes, w) array with x,y coordinates
        self.edges = None  # (n_edges, 2) array with node indices
        self.triangles = None  # (n_triangles, 3) array with node indices
        self.domains = []  # List of Domain objects

    def read_from_file(self, filename: str):
        with open(filename, 'r') as f:
            data = f.read()

        lines = data.strip().split('\n')

        # Parse nodes
        i = 0
        line = lines[i]
        n_nodes = int(line.split()[-1])
        self.nodes = np.zeros((n_nodes, 2))

        for j in range(n_nodes):
            i += 1
            if i >= len(lines) or "Number of" in lines[i]:
                break
            parts = lines[i].split(':')
            if len(parts) < 2:
                continue

            node_id = int(parts[0])
            coords = parts[1].strip().split()
            x = float(coords[0])
            y = float(coords[1])
            self.nodes[node_id] = [x, y]

        # Find and parse edges
        while i < len(lines):
            if "Number of edges" in lines[i]:
                n_edges = int(lines[i].split()[-1])
                self.edges = np.zeros((n_edges, 2), dtype=int)

                for j in range(n_edges):
                    i += 1
                    if i >= len(lines) or "Number of" in lines[i]:
                        break
                    parts = lines[i].split(':')
                    if len(parts) < 2:
                        continue

                    edge_id = int(parts[0])
                    nodes = parts[1].strip().split()
                    n1 = int(nodes[0])
                    n2 = int(nodes[1])
                    self.edges[edge_id] = [n1, n2]

            # Parse triangles
            elif "Number of triangles" in lines[i]:
                n_triangles = int(lines[i].split()[-1])
                self.triangles = np.zeros((n_triangles, 3), dtype=int)

                for j in range(n_triangles):
                    i += 1
                    if i >= len(lines) or "Number of" in lines[i]:
                        break
                    parts = lines[i].split(':')
                    if len(parts) < 2:
                        continue

                    tri_id = int(parts[0])
                    nodes = parts[1].strip().split()
                    n1 = int(nodes[0])
                    n2 = int(nodes[1])
                    n3 = int(nodes[2])
                    self.triangles[tri_id] = [n1, n2, n3]

            self.domains = []
            i += 1

    def write_to_file(self, filename: str):
        with open(filename, 'w') as f:
            # Write nodes
            f.write(f"Number of nodes {len(self.nodes)}\n")
            for i, node in enumerate(self.nodes):
                f.write(f"{i:6d} : {node[0]:15.7e} {node[1]:15.7e}\n")

            # Write edges
            f.write(f"Number of edges {len(self.edges)}\n")
            for i, edge in enumerate(self.edges):
                f.write(f"{i:6d} : {edge[0]:6d} {edge[1]:6d}\n")

            # Write triangles
            if self.triangles is not None:
                f.write(f"Number of triangles {len(self.triangles)}\n")
                for i, tri in enumerate(self.triangles):
                    f.write(f"{i:6d} : {tri[0]:6d} {tri[1]:6d} {tri[2]:6d}\n")

            # Write domains
            f.write(f"Number of domains {len(self.domains)}\n")
            for i, domain in enumerate(self.domains):
                f.write(f"  Domain : {i:6d}\n")
                f.write(f"  Name : {domain.name}\n")
                f.write(f"  Number of elements : {len(domain.elements):6d}\n")

                # Write elements with 10 per line
                elements_per_line = 10
                for j in range(0, len(domain.elements), elements_per_line):
                    end = min(j + elements_per_line, len(domain.elements))
                    # elements_str = "    " + " ".join(f"{elem:6d}" for elem in domain.elements[j:end])
                    elements_str = " ".join(f"{elem:6d}" for elem in domain.elements[j:end])
                    f.write(f"{elements_str}\n")

    def find_pillars(self):
        """Identify bridge pillars based on vertical edge structures"""
        # Find the minimum y-coordinate (bottom of bridge)
        min_y = np.min(self.nodes[:, 1])

        # Find all horizontal edges at the bottom of the structure
        bottom_edges = []
        for i, (n1, n2) in enumerate(self.edges):
            y1, y2 = self.nodes[n1, 1], self.nodes[n2, 1]

            # If both nodes are at the bottom and edge is horizontal
            if abs(y1 - min_y) < 1e-5 and abs(y2 - min_y) < 1e-5:
                bottom_edges.append(i)

        # Group edges into pillars by connectivity
        pillars = []
        visited = set()

        for edge_idx in bottom_edges:
            if edge_idx in visited:
                continue

            # Start a new pillar
            current_pillar = []
            edge_stack = [edge_idx]

            while edge_stack:
                current_edge = edge_stack.pop()
                if current_edge in visited:
                    continue

                visited.add(current_edge)
                current_pillar.append(current_edge)

                # Get nodes of this edge
                n1, n2 = self.edges[current_edge]

                # Find connected edges at the bottom
                for i, (m1, m2) in enumerate(self.edges):
                    if i in visited or i not in bottom_edges:
                        continue

                    # If edges share a node and both at bottom
                    if n1 == m1 or n1 == m2 or n2 == m1 or n2 == m2:
                        edge_stack.append(i)

            # If we found a pillar, add it
            if current_pillar:
                pillars.append(current_pillar)

        return pillars

    def get_pillar_sides(self, pillar_bottom_edges, side_height_factor=0.5):
        """Find side edges for a given pillar up to a specified height"""
        # Get all nodes in the pillar bottom
        pillar_nodes = set()
        for edge_idx in pillar_bottom_edges:
            n1, n2 = self.edges[edge_idx]
            pillar_nodes.add(n1)
            pillar_nodes.add(n2)

        # Calculate pillar width
        x_values = [self.nodes[n, 0] for n in pillar_nodes]
        pillar_width = max(x_values) - min(x_values)

        # Calculate maximum height for side edges
        min_y = min(self.nodes[n, 1] for n in pillar_nodes)
        max_height = min_y + side_height_factor * pillar_width

        # Find vertical or near-vertical edges connected to the pillar bottom
        side_edges = []

        for i, (n1, n2) in enumerate(self.edges):
            # Skip if already in pillar bottom
            if i in pillar_bottom_edges:
                continue

            # # Check if edge is connected to pillar bottom
            # if n1 not in pillar_nodes and n2 not in pillar_nodes:
            #     continue
            if self.nodes[n1, 0] > max(x_values) * 1.01 or self.nodes[n2, 0] < min(x_values) * 0.99:
                continue

            y1, y2 = self.nodes[n1, 1], self.nodes[n2, 1]

            # Include edge if it's within the height limit
            if min(y1, y2) <= max_height and max(y1, y2) <= max_height:
                side_edges.append(i)

        return side_edges

    def find_extremity_corners(self):
        """Find lower corners of left and right extremities"""
        # Find min/max x coordinates
        min_x = np.min(self.nodes[:, 0])
        max_x = np.max(self.nodes[:, 0])
        max_y = np.max(self.nodes[:, 1])

        left_corner_edges = []
        right_corner_edges = []

        min_corn_y = 1e12

        for i, (n1, n2) in enumerate(self.edges):
            x1, y1 = self.nodes[n1, 0], self.nodes[n1, 1]
            x2, y2 = self.nodes[n2, 0], self.nodes[n2, 1]

            # Left corner
            if abs(x1 - min_x) < 1e-5 and abs(x2 - min_x) < 1e-5:
                left_corner_edges.append(i)
                min_corn_y = min(min_corn_y, min(y1, y2))

            # Right corner
            if abs(x1 - max_x) < 1e-5 and abs(x2 - max_x) < 1e-5:
                right_corner_edges.append(i)
                min_corn_y = min(min_corn_y, min(y1, y2))

        extremity_height = (max_y - min_corn_y)


        for i, (n1, n2) in enumerate(self.edges):
            x1, y1 = self.nodes[n1, 0], self.nodes[n1, 1]
            x2, y2 = self.nodes[n2, 0], self.nodes[n2, 1]

            # Left corner
            if abs(x1 - min_x) < extremity_height and abs(x2 - min_x) < extremity_height:
                if abs(y1 - min_corn_y) < 1e-5 and abs(y2 - min_corn_y) < 1e-5:
                    left_corner_edges.append(i)

            # Right corner
            if abs(x1 - max_x) < extremity_height and abs(x2 - max_x) < extremity_height:
                if abs(y1 - min_corn_y) < 1e-5 and abs(y2 - min_corn_y) < 1e-5:
                    right_corner_edges.append(i)



        return left_corner_edges, right_corner_edges

    def add_boundary_domains(self):
        """Add domains for the boundary conditions as described"""
        # Find pillars
        pillars = self.find_pillars()

        # Create domains for pillar bottoms
        next_domain_id = len(self.domains)

        # Add pillar domains
        for i, pillar_bottom in enumerate(pillars):
            # Get side edges up to half the pillar width
            side_edges = self.get_pillar_sides(pillar_bottom)

            # Combine bottom and side edges
            all_edges = pillar_bottom + side_edges

            # Create domain
            domain_name = f"PillarBottom_{i}"
            self.domains.append(Domain(domain_name, np.array(all_edges, dtype=int)))

        # Add extremity corners
        left_corner, right_corner = self.find_extremity_corners()
        print(f"Left corner edges: {left_corner}")
        print(f"Right corner edges: {right_corner}")

        if left_corner:
            self.domains.append(Domain("LeftCorner", np.array(left_corner, dtype=int)))

        if right_corner:
            self.domains.append(Domain("RightCorner", np.array(right_corner, dtype=int)))

        return next_domain_id, next_domain_id + len(pillars) + 2

    def visualize_domain(self, domain_idx, fig, ax, title=None):
        """Visualize a specific domain"""

        # Draw all edges lightly
        lines = []
        for n1, n2 in self.edges:
            lines.append([self.nodes[n1, :2], self.nodes[n2, :2]])

        lc = LineCollection(lines, colors='lightgray', linewidths=0.5)
        ax.add_collection(lc)

        # Highlight domain edges
        if domain_idx < len(self.domains):
            domain = self.domains[domain_idx]
            domain_lines = []

            for edge_idx in domain.elements:
                if edge_idx < len(self.edges):
                    n1, n2 = self.edges[edge_idx]
                    domain_lines.append([self.nodes[n1, :2], self.nodes[n2, :2]])

            lc2 = LineCollection(domain_lines, colors='red', linewidths=2)
            ax.add_collection(lc2)

            # if title:
                # ax.set_title(f"{title}: {domain.name}")
            # else:
                # ax.set_title(f"Domain {domain_idx}: {domain.name}")

        # ax.set_aspect('equal')
        # plt.autoscale()
        # plt.tight_layout()
        # plt.show()

def main():
    if len(sys.argv) < 3:
        print("Usage: python add_bridge_boundary_domains.py input_mesh.txt output_mesh.txt")
        return

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Load mesh
    mesh = Mesh()
    mesh.read_from_file(input_file)

    # Add boundary domains
    start_idx, end_idx = mesh.add_boundary_domains()

    print(f"Added {end_idx - start_idx} new domains for boundary conditions")

    # Visualize the domains
    fig, ax = plt.subplots(figsize=(10, 8))
    for i in range(start_idx, end_idx):
        if i < len(mesh.domains):
            mesh.visualize_domain(i, fig, ax, "New Boundary Domain")
    plt.show()

    # Save the updated mesh
    mesh.write_to_file(output_file)
    print(f"Updated mesh saved to {output_file}")

if __name__ == "__main__":
    main()
