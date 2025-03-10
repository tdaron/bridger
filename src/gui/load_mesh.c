#include "../mesh/mesh.h"
#include "gui.h"

/**
@returns This functions returns the VBO of the mesh.
**/
unsigned int load_mesh_into_vao(unsigned int VAO, Mesh *mesh) {
  glBindVertexArray(VAO);

  unsigned int VBO;
  glGenBuffers(1, &VBO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);

  // Loading mesh vertices into the vbo
  glBufferData(GL_ARRAY_BUFFER, mesh->vertexArraySize * sizeof(float),
               mesh->vertexArray, GL_STATIC_DRAW);

  // 3 coordinates per vertex
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
  glEnableVertexAttribArray(0);

  return VBO;
}

/**
@returns This functions returns the VBO of the mesh.
**/
unsigned int load_field_into_vao(unsigned int VAO, float *field, int size) {
  glBindVertexArray(VAO);

  unsigned int VBO;
  glGenBuffers(1, &VBO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);

  // Loading mesh vertices into the vbo
  glBufferData(GL_ARRAY_BUFFER, size * sizeof(float), field, GL_STATIC_DRAW);

  // 3 coordinates per vertex
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
  glEnableVertexAttribArray(1);

  return VBO;
}
