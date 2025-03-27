#include "gui/gui.h"
#include "mesh/mesh.h"
#include <gmshc.h>
#include <log.h>
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
double scale = 0;
unsigned int meshVBO;
unsigned int fieldVBO;
int windowWidth, windowHeight;
float tankDx = 0;
int main() {
  int ierr;
  ss_init();
  gmshInitialize(0, NULL, 1, 0, &ierr);
  printf("GMSH loaded\n");

  mesh = generate_mesh(&settings, &scale);
  field = compute_field(mesh, &settings);
  // Mesh *mesh = readMesh("data/mesh.txt");

  GLFWwindow *window = createWindow(800, 600, "Bridger");
  unsigned int shaderProgram = get_program("src/gui/shaders/vertex.vert",
                                           "src/gui/shaders/fragment.frag");

  vao = create_vao();
  meshVBO = load_mesh_into_vao(vao, mesh, 0);
  fieldVBO = load_field_into_vao(vao, field, mesh->numTriangles * 3 * 3, 0);

  glfwSetMouseButtonCallback(window, mouse_callback);

  time_t old = time(NULL);
  while (!glfwWindowShouldClose(window)) {
    time_t new = time(NULL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    processInput(window);

    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glUseProgram(shaderProgram);
    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLES, 0, mesh->numTriangles * 3);
    DrawImage("assets/tank.png", -0.9+tankDx, 0.2, 0.8);
    glfwGetWindowSize(window, &windowWidth, &windowHeight);
    glfwSwapBuffers(window);
    glfwPollEvents();
    // glfwGetCursorPos(window, &xpos, &ypos);
    // printf("%f %f\n", xpos, ypos );
  }
  glfwTerminate();
  freeMesh(mesh);
  return 0;
}

void processInput(GLFWwindow *window) {
  if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    glfwSetWindowShouldClose(window, true);
  if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) tankDx += 0.01;
  if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) tankDx -= 0.01;
}

void get_cursor_position(GLFWwindow *window, MeshSettings *s, double *xpos,
                         double *ypos) {

  glfwGetCursorPos(window, xpos, ypos);
  log_info("Click on: %f %f - Width: %d Height: %d", *xpos, *ypos, windowWidth,
           windowHeight);

  // Normalize to range [-1, 1]
  *xpos = 2.0 * (*xpos / windowWidth) - 1.0;
  *ypos = 1.0 - 2.0 * (*ypos / windowHeight);
  log_info("After normalization: %f %f", *xpos, *ypos);
}

void mouse_callback(GLFWwindow *window, int button, int action, int mods) {
  // Left button  &&  pushed
  if (button == 0 && action == 0) {
    double xpos;
    double ypos;
    get_cursor_position(window, &settings, &xpos, &ypos);
    settings.holeX = dn(xpos, mesh, 0);
    settings.holeY = dn(ypos, mesh, 1);
    freeMesh(mesh);
    mesh = generate_mesh(&settings, &scale);
    printf("Adding hole at %f %f\n", dn(xpos, mesh, 0), dn(ypos, mesh, 1));
    free(field);
    field = compute_field(mesh, &settings);
    meshVBO = load_mesh_into_vao(vao, mesh, meshVBO);
    fieldVBO =
        load_field_into_vao(vao, field, mesh->numTriangles * 3 * 3, fieldVBO);
    glfwPostEmptyEvent();
  }
}
