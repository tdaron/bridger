import numpy as np
import matplotlib.pyplot as plt

def compute_norms(filename):
    norms = []
    with open(filename, 'r') as file:
        for line in file:
            values = line.split()
            if len(values) == 2:  # Ensure there are exactly two columns
                try:
                    x, y = map(float, values)
                    norm = np.sqrt(x**2 + y**2)
                    norms.append(norm)
                except ValueError:
                    pass  # Ignore lines that cannot be converted
    return norms

def plot_histogram(norms, bins=200):
    plt.hist(norms, bins=bins, edgecolor='black', alpha=0.7)
    plt.xlabel('Norm Values')
    plt.ylabel('Frequency')
    plt.title('Histogram of Norms')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()

# Example usage
filename = './data/solution.txt'  # Change this to the actual file path
norms = compute_norms(filename)
plot_histogram(norms)
