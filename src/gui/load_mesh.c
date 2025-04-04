#include "../mesh/mesh.h"
#include "gui.h"
#include "problem.h"
#include <log.h>

/**
@returns This functions returns the VBO of the mesh.
**/
unsigned int load_mesh_into_vao(unsigned int VAO, Mesh *mesh,
                                unsigned int VBO) {
  glBindVertexArray(VAO);

  if (VBO == 0) {
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);

    // Loading mesh vertices into the vbo
    glBufferData(GL_ARRAY_BUFFER, mesh->vertexArraySize * sizeof(float),
                 mesh->vertexArray, GL_STATIC_DRAW);
  } else {
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0, mesh->vertexArraySize * sizeof(float),
                    mesh->vertexArray);
  }
  // 3 coordinates per vertex
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
  glEnableVertexAttribArray(0);

  return VBO;
}

/**
@returns This functions returns the VBO of the mesh.
**/
unsigned int load_field_into_vao(unsigned int VAO, double *field, int size,
                                 unsigned int VBO) {
  glBindVertexArray(VAO);

  if (VBO == 0) {
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    // Loading mesh vertices into the vbo
    glBufferData(GL_ARRAY_BUFFER, size * sizeof(double), field, GL_STATIC_DRAW);

  } else {
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0, size * sizeof(double), field);
  }

  // 3 coordinates per vertex
  glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double),
                        (void *)0);
  glEnableVertexAttribArray(1);

  return VBO;
}

/**
@returns This functions returns the VBO of the mesh.
**/
unsigned int load_soluce_into_vao(unsigned int VAO, double *soluce,
                                  geo *geometry, unsigned int VBO,
                                  double scale) {
  glBindVertexArray(VAO);

  size_t size = geometry->theElements->nElem * 3 *
                2; // each triangle is 3 nodes 2 dimensions
  double *field = malloc(sizeof(double) * size);
  for (size_t i = 0; i < geometry->theElements->nElem; i++) {
    // For each triangle
    for (size_t j = 0; j < 3; j++) {
      // 3 vertices
      int node_id = geometry->theElements->elem[i * 3 + j];

      for (size_t k = 0; k < 2; k++) {
        field[i * 3 * 2 + j * 2 + k] = soluce[node_id * 2 + k];
      }
    }
  }

  if (VBO == 0) {
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    // Loading mesh vertices into the vbo
    glBufferData(GL_ARRAY_BUFFER, size * sizeof(double), field, GL_STATIC_DRAW);

  } else {
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0, size * sizeof(double), field);
  }

  // 2 coordinates per vertex
  glVertexAttribPointer(2, 2, GL_DOUBLE, GL_FALSE, 2 * sizeof(double),
                        (void *)0);
  glEnableVertexAttribArray(2);

  free(field);
  free(soluce);
  return VBO;
}
