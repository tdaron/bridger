#include "gui/gui.h"
#include "mesh/mesh.h"
#include <stdbool.h>
#include <gmshc.h>
#include <stdio.h>

void processInput(GLFWwindow *window);

int main() {
  int ierr;
  gmshInitialize(0, NULL, 1, 0, &ierr);
  printf("GMSH loaded\n");
  GLFWwindow *window = createWindow(800, 600, "Bridger");

  Mesh *mesh = readMesh("data/mesh.txt");

  unsigned int shaderProgram = get_program("shaders/vertex.vert", "shaders/fragment.frag");


  unsigned int vao = create_vao();
  unsigned int meshVBO = load_mesh_into_vao(vao, mesh);

  while (!glfwWindowShouldClose(window)) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    processInput(window);

    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glUseProgram(shaderProgram);
    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLES, 0, mesh->numTriangles * 3);

    glfwSwapBuffers(window);
    glfwPollEvents();
  }
  glfwTerminate();
  freeMesh(mesh);
  return 0;
}

void processInput(GLFWwindow *window) {
  if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    glfwSetWindowShouldClose(window, true);
}
