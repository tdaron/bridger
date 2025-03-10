#include "gui/gui.h"
#include "mesh/mesh.h"
#include <gmshc.h>
#include <solver.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

void processInput(GLFWwindow *window);
void get_cursor_position(GLFWwindow *window, MeshSettings *s, double *xpos,
                         double *ypos);
void mouse_callback(GLFWwindow *window, int button, int action, int mods);

MeshSettings settings = {
    .bridgeHeight = .8,
    .bridgeWidth = 20,
    .pillarsWidth = 1.5,
    .pillarsHeight = 3.5,
    .pillarsNumber = 4,
    .offset = 3,
    .baseElementSize = 0.5,
    .preciseElementSize = 0.1,
    .precisionRadius = 2,
};

Mesh *mesh;
unsigned int vao;
float *field;
int main() {
  int ierr;
  ss_init();
  gmshInitialize(0, NULL, 1, 0, &ierr);
  printf("GMSH loaded\n");

  mesh = generate_mesh(&settings);
  field = compute_field(mesh, &settings);
  // Mesh *mesh = readMesh("data/mesh.txt");

  GLFWwindow *window = createWindow(800, 600, "Bridger");
  unsigned int shaderProgram = get_program("src/gui/shaders/vertex.vert",
                                           "src/gui/shaders/fragment.frag");

  vao = create_vao();
  unsigned int meshVBO = load_mesh_into_vao(vao, mesh);
  unsigned int fieldVBO =
      load_field_into_vao(vao, field, mesh->numTriangles * 3 * 3);

  glfwSetMouseButtonCallback(window, mouse_callback);

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

void get_cursor_position(GLFWwindow *window, MeshSettings *s, double *xpos,
                         double *ypos) {
  int width, height;
  glfwGetFramebufferSize(window, &width, &height);

  double dw = (double)width;
  double dh = (double)height;

  glfwGetCursorPos(window, xpos, ypos);

  // Normalize to range [-1, 1]
  *xpos = 2.0 * (*xpos / dw) - 1.0;
  *ypos = 2.0 * (*ypos / dh) - 1.0;
  *ypos *= -1;
}

void mouse_callback(GLFWwindow *window, int button, int action, int mods) {
  // Left button  &&  pushed
  if (button == 0 && action == 0) {
    double xpos;
    double ypos;
    get_cursor_position(window, &settings, &xpos, &ypos);
    settings.holeX = dn(xpos, mesh, 0);
    settings.holeY = dn(ypos, mesh, 1);
    mesh = generate_mesh(&settings);
    printf("Adding hole at %f %f\n", dn(xpos, mesh, 0), dn(ypos, mesh, 1));
    field = compute_field(mesh, &settings);
    unsigned int meshVBO = load_mesh_into_vao(vao, mesh);
    unsigned int fieldVBO =
        load_field_into_vao(vao, field, mesh->numTriangles * 3 * 3);
  }
}
