import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from typing import List, Dict, Tuple
import random
from datetime import datetime

class Domain:
    def __init__(self, name: str, elements: np.ndarray):
        self.name = name
        self.elements = elements

    def __repr__(self):
        return f"Domain(name='{self.name}', elements={self.elements})"

class Mesh:
    def __init__(self):
        self.nodes = None  # (n_nodes, 2) array with x,y coordinates
        self.edges = None  # (n_edges, 2) array with node indices
        self.elements = {}  # Dictionary of element types: {'triangles': np.array, 'quads': np.array}
        self.domains = []  # List of Domain objects

    def read_from_string(self, data: str):
        lines = data.strip().split('\n')

        # Parse nodes
        i = 0
        line = lines[i]
        n_nodes = int(line.split()[-1])
        self.nodes = np.zeros((n_nodes, 2))

        for j in range(n_nodes):
            i += 1
            if i >= len(lines) or "Number of" in lines[i]:
                break  # We've reached the end of nodes section
            parts = lines[i].split()
            if len(parts) < 4:  # Skip lines that don't match the expected format
                continue
            node_id = int(parts[0])
            x = float(parts[2])
            y = float(parts[3])
            self.nodes[node_id] = [x, y]

        # Find and parse edges
        while i < len(lines):
            if "Number of edges" in lines[i]:
                n_edges = int(lines[i].split()[-1])
                self.edges = np.zeros((n_edges, 2), dtype=int)

                for j in range(n_edges):
                    i += 1
                    if i >= len(lines) or "Number of" in lines[i]:
                        break  # We've reached the end of edges section
                    parts = lines[i].split()
                    if len(parts) < 4:  # Skip lines that don't match the expected format
                        continue
                    edge_id = int(parts[0])
                    n1 = int(parts[2])
                    n2 = int(parts[3])
                    self.edges[edge_id] = [n1, n2]

            # Parse triangles
            elif "Number of triangles" in lines[i]:
                n_triangles = int(lines[i].split()[-1])
                triangles = np.zeros((n_triangles, 3), dtype=int)

                for j in range(n_triangles):
                    i += 1
                    if i >= len(lines) or "Number of" in lines[i]:
                        break  # We've reached the end of triangles section
                    parts = lines[i].split()
                    if len(parts) < 5:  # Skip lines that don't match the expected format
                        continue
                    tri_id = int(parts[0])
                    n1 = int(parts[2])
                    n2 = int(parts[3])
                    n3 = int(parts[4])
                    triangles[tri_id] = [n1, n2, n3]

                self.elements['triangles'] = triangles

            # Parse quads
            elif "Number of quads" in lines[i]:
                n_quads = int(lines[i].split()[-1])
                quads = np.zeros((n_quads, 4), dtype=int)

                for j in range(n_quads):
                    i += 1
                    if i >= len(lines) or "Number of" in lines[i]:
                        break  # We've reached the end of quads section
                    parts = lines[i].split()
                    if len(parts) < 6:  # Skip lines that don't match the expected format
                        continue
                    quad_id = int(parts[0])
                    n1 = int(parts[2])
                    n2 = int(parts[3])
                    n3 = int(parts[4])
                    n4 = int(parts[5])
                    quads[quad_id] = [n1, n2, n3, n4]

                self.elements['quads'] = quads

            # Parse domains
            elif "Number of domains" in lines[i]:
                n_domains = int(lines[i].split()[-1])

                d = 0
                while d < n_domains and i < len(lines):
                    i += 1
                    if "Domain :" in lines[i]:
                        domain_id = int(lines[i].split()[-1])

                        i += 1
                        domain_name = lines[i].split(':', 1)[1].strip()

                        i += 1
                        n_elements = int(lines[i].split()[-1])

                        # Read elements (they may span multiple lines)
                        elements = []
                        while len(elements) < n_elements and i + 1 < len(lines):
                            i += 1
                            if "Domain :" in lines[i]:  # We've hit the next domain
                                i -= 1  # Go back one line
                                break
                            elements.extend(int(x) for x in lines[i].split())

                        domain_elements = np.array(elements[:n_elements], dtype=int)
                        self.domains.append(Domain(domain_name, domain_elements))
                        d += 1
                break  # We've parsed all domains, exit the while loop

            i += 1

    def get_domain_centers(self):
        """Calculate the center position for each domain's label"""
        domain_centers = {}

        for domain in self.domains:
            # Get all nodes that make up the edges in this domain
            all_nodes = []
            for edge_idx in domain.elements:
                if edge_idx < len(self.edges):
                    node_indices = self.edges[edge_idx]
                    all_nodes.extend(node_indices)

            # Get unique nodes
            unique_nodes = np.unique(all_nodes)

            # Calculate the center of these nodes
            if len(unique_nodes) > 0:
                center_x = np.mean(self.nodes[unique_nodes, 0])
                center_y = np.mean(self.nodes[unique_nodes, 1])
                domain_centers[domain.name] = (center_x, center_y)

        return domain_centers

    def visualize(self, current_time=None, user=None):
        fig, ax = plt.subplots(figsize=(12, 10))

        # Draw all elements as background with light gray
        for element_type, elements in self.elements.items():
            for element in elements:
                # Extract coordinates for each node in the element
                element_coords = self.nodes[element]

                # Plot the element as a light gray polygon
                poly = plt.Polygon(element_coords, alpha=0.2, color='lightgray',
                                  edgecolor='gray', linewidth=0.5)
                ax.add_patch(poly)

        # Get domain colors - use a colormap for distinct colors
        colors = plt.cm.tab10(np.linspace(0, 1, len(self.domains)))
        random.shuffle(colors)  # Mix colors for better contrast

        # Create a dictionary to store edges by domain for the legend
        domain_lines = {}
        domain_colors = {}

        # Plot domains (each domain has edges with different colors)
        for i, domain in enumerate(self.domains):
            # Get the edge indices for this domain
            edge_indices = domain.elements

            # Lines for this domain
            lines = []
            for edge_idx in edge_indices:
                if edge_idx >= len(self.edges):
                    continue  # Skip if edge index is out of bounds

                # Get the nodes of this edge
                node_indices = self.edges[edge_idx]

                # Extract coordinates
                line_coords = self.nodes[node_indices]
                lines.append(line_coords)

            # Create line collection for this domain
            line_collection = LineCollection(lines, color=colors[i], linewidth=2.5, label=domain.name)
            domain_lines[domain.name] = (line_collection, len(edge_indices), colors[i])
            domain_colors[domain.name] = colors[i]
            ax.add_collection(line_collection)

        # Calculate domain centers for labels
        domain_centers = self.get_domain_centers()

        # Add domain labels directly on the plot
        for domain_name, (center_x, center_y) in domain_centers.items():
            color = domain_colors.get(domain_name, 'black')

            # Add label with a white background for readability
            ax.text(center_x, center_y, domain_name,
                    fontsize=9, weight='bold', color=color,
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor=color, boxstyle='round,pad=0.3'),
                    ha='center', va='center', zorder=10)

        # Plot all nodes
        ax.scatter(self.nodes[:, 0], self.nodes[:, 1], c='black', s=5, zorder=5)

        # Plot domain legends
        handles = []
        labels = []
        for name, (_, count, color) in sorted(domain_lines.items()):
            line = plt.Line2D([0], [0], color=color, linewidth=2.5)
            handles.append(line)
            labels.append(f"{name} ({count} edges)")

        ax.legend(handles, labels, loc='upper right')

        # Set axis properties
        ax.set_aspect('equal')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')

        # Create title based on element types
        element_types = []
        for elem_type, elems in self.elements.items():
            element_types.append(f"{len(elems)} {elem_type}")

        element_info = ", ".join(element_types)
        ax.set_title(f'Mesh Visualization with {len(self.domains)} Domains\n({element_info})', fontsize=14)

        # Add metadata at bottom
        if current_time is None:
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if user is None:
            user = "User"

        plt.figtext(0.01, 0.01, f"Created: {current_time} by {user}", fontsize=8)

        # Set axis limits from data with a small margin
        x_min, x_max = np.min(self.nodes[:, 0]), np.max(self.nodes[:, 0])
        y_min, y_max = np.min(self.nodes[:, 1]), np.max(self.nodes[:, 1])
        margin = 0.05 * max(x_max - x_min, y_max - y_min)
        ax.set_xlim(x_min - margin, x_max + margin)
        ax.set_ylim(y_min - margin, y_max + margin)

        plt.tight_layout()
        plt.savefig('mesh_visualization.png', dpi=300, bbox_inches='tight')
        plt.show()

        print(f"Visualization saved as 'mesh_visualization.png'")

def read_mesh_from_file(filename):
    with open(filename, 'r') as f:
        data = f.read()

    mesh = Mesh()
    mesh.read_from_string(data)
    return mesh

def main():
    # Import from file
    with open('data/mesh10.txt', 'r') as f:
        mesh_data = f.read()

    # Create a mesh object
    mesh = Mesh()

    # Parse the mesh data
    mesh.read_from_string(mesh_data)

    # Visualize the mesh with current time and user
    current_time = "2025-03-25 20:46:47"
    user = "mxdbck"
    mesh.visualize(current_time, user)

if __name__ == "__main__":
    main()
