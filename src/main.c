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

void regen_mesh();
void regen_solution();
MeshSettings settings = {.bridgeHeight = .8,
                         .bridgeWidth = 20,
                         .pillarsWidth = 1.5,
                         .pillarsHeight = 3.5,
                         .pillarsNumber = 4,
                         .offset = 3,
                         .baseElementSize = 0.5,
                         .preciseElementSize = 0.1,
                         .precisionRadius = 2,
                         .tankLength = 4,
                         .tankWeight = -100000,
                         .tankX = 0};

Mesh *gpu_mesh;
unsigned int vao;
double *field;
double scale = 0;
unsigned int meshVBO;
unsigned int fieldVBO;
unsigned int soluceVBO;
int windowWidth, windowHeight;
float tankDx = 0;
geo *geometry;

int seeDisformation = 1.0;

int main() {
  int ierr;
  ss_init();
  gmshInitialize(0, NULL, 1, 0, &ierr);
  printf("GMSH loaded\n");

  GLFWwindow *window = createWindow(800, 600, "Bridger");
  unsigned int shaderProgram = get_program("src/gui/shaders/vertex.vert",
                                           "src/gui/shaders/fragment.frag");

  vao = create_vao();
  regen_mesh();
  settings.tankX = ((-0.9 + tankDx) / scale) + gpu_mesh->centers[0];
  regen_solution();

  glfwSetMouseButtonCallback(window, mouse_callback);

  time_t old = time(NULL);

  unsigned int seeDisformationUniform =
      glGetUniformLocation(shaderProgram, "seeDisformation");

  while (!glfwWindowShouldClose(window)) {
    time_t new = time(NULL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    processInput(window);

    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glUseProgram(shaderProgram);
    glUniform1i(seeDisformationUniform, seeDisformation);
    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLES, 0, gpu_mesh->numTriangles * 3);
    DrawImage("assets/winki.png", -0.9 + tankDx, 0.2, 0.8);
    glfwGetWindowSize(window, &windowWidth, &windowHeight);
    glfwSwapBuffers(window);
    settings.tankX = ((-0.9 + tankDx) / scale) + gpu_mesh->centers[0];
    glfwPollEvents();
    // glfwGetCursorPos(window, &xpos, &ypos);
    // printf("%f %f\n", xpos, ypos );
  }
  glfwTerminate();
  freeMesh(gpu_mesh);
  return 0;
}

static int prevTKeyState = GLFW_RELEASE; // Store previous state of T key
void processInput(GLFWwindow *window) {
  if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    glfwSetWindowShouldClose(window, true);
  if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) {
    tankDx += 0.01;
    regen_solution();
  }
  if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) {
    tankDx -= 0.01;
    regen_solution();
  }
  if (glfwGetKey(window, GLFW_KEY_T) == GLFW_PRESS &&
      prevTKeyState == GLFW_RELEASE) {
    seeDisformation = !seeDisformation; // Toggle only on key press (not hold)
  }

  prevTKeyState = glfwGetKey(window, GLFW_KEY_T); // Update previous state
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

void regen_mesh() {
  freeMesh(gpu_mesh);
  gpu_mesh = generate_mesh(&settings, &scale);
  free(field);
  field = compute_field(gpu_mesh, &settings);
  meshVBO = load_mesh_into_vao(vao, gpu_mesh, meshVBO);
  fieldVBO =
      load_field_into_vao(vao, field, gpu_mesh->numTriangles * 3 * 3, fieldVBO);
  regen_solution();
}

void regen_solution() {

  int *tankEdges;
  int nTankEdges;
  find_tank_edges(gpu_mesh, &tankEdges, &nTankEdges, &settings);
  double *soluce = compute_solution("solver/data/mesh.txt", NULL, &geometry, nTankEdges,
                                    tankEdges, settings.tankWeight);
  soluceVBO = load_soluce_into_vao(vao, soluce, geometry, 0, gpu_mesh->scale);
}

void mouse_callback(GLFWwindow *window, int button, int action, int mods) {
  // Left button  &&  pushed
  if (button == 0 && action == 0) {
    // double xpos;
    // double ypos;
    // get_cursor_position(window, &settings, &xpos, &ypos);
    // settings.holeX = dn(xpos, gpu_mesh, 0);
    // settings.holeY = dn(ypos, gpu_mesh, 1);
    // printf("Adding hole at %f %f\n", dn(xpos, gpu_mesh, 0),
    //        dn(ypos, gpu_mesh, 1));
    regen_mesh();
    glfwPostEmptyEvent();
  }
}
