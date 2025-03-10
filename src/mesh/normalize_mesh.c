#include "mesh.h"
#include <log.h>

void normalize_mesh(Mesh *mesh) {
  log_info("Normalizing mesh coordinates: num triangles: %d", mesh->numTriangles);
  float max[3] = {0, 0, 0};
  float min[3] = {0, 0, 0};
  for (int i = 0; i < mesh->numTriangles; i++) {
    // For each triangle
    for (int j = 0; j < 3; j++) {
      // For each vertex
      for (int k = 0; k < 3; k++) {
        float value = mesh->vertexArray[i * 9 + j * 3 + k];
        if (value > max[k])
          max[k] = value;
        if (value < min[k])
          min[k] = value;
      }
    }
  }

  // width, height, depth
  float deltas[3] = {max[0] - min[0], max[1] - min[1], max[2] - min[2]};

  float scale =
      1.9f / (deltas[0] > deltas[1]
                  ? deltas[0]
                  : deltas[1]); // Scale width or height to fit in [-0.95, 0.95]

  float centers[3] = {(min[0] + max[0]) / 2, (min[1] + max[1]) / 2,
                      (min[2] + max[2]) / 2};
  log_info("Mesh bounds: X[%.2f, %.2f, %.2f], Y[%.2f, %.2f, %.2f]", min[0],
           max[0], min[1], max[1], min[2], max[2]);
  log_info("Center: (%.2f, %.2f, %.2f), Scale: %.5f", centers[0], centers[1],
           centers[2], scale);

  mesh->centers[0] = centers[0];
  mesh->centers[1] = centers[1];
  mesh->centers[2] = centers[2];

  mesh->scale = scale;
  for (int i = 0; i < mesh->numTriangles; i++) {
    // For each triangle
    for (int j = 0; j < 3; j++) {
      // For each vertex
      for (int k = 0; k < 3; k++) {
        #define val(i, j, k) mesh->vertexArray[i * 9 + j * 3 + k]
        val(i, j, k) = (val(i, j, k) - centers[k]) * scale;
      }
    }
  }
}

double dn(double coord, Mesh* mesh, int k) {
  return (coord / mesh->scale) + mesh->centers[k];
}
