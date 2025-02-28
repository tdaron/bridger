// clang-format off
#include <glad/glad.h>
#include <GLFW/glfw3.h>
// clang-format on
#include "gui/gui.h"
#include "mesh/read_mesh.h"
#include <stdbool.h>

void processInput(GLFWwindow *window);

int main() {
  GLFWwindow *window = createWindow(800, 600, "Bridger");

  float vertices[] = {
      // x y z r g b
      -0.5f, -0.5f, 0.0f, 1.0f, 0, 0, // first triangle
      0.5f,  -0.5f, 0.0f, 0,    1, 0, // second triangle
      0.0f,  0.5f,  0.0f, 0,    0, 1  // third triangle
  };

  unsigned int VBO;
  glGenBuffers(1, &VBO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  Mesh *mesh = readMeshForOpenGL("data/mesh.txt");

  // glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

  unsigned int vertexShader =
      compile_shader("shaders/vertex.vert", GL_VERTEX_SHADER);
  unsigned int fragmentShader =
      compile_shader("shaders/fragment.frag", GL_FRAGMENT_SHADER);
  unsigned int shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, vertexShader);
  glAttachShader(shaderProgram, fragmentShader);
  glLinkProgram(shaderProgram);
  glDeleteShader(vertexShader);
  glDeleteShader(fragmentShader);

  unsigned int VAO;
  glGenVertexArrays(1, &VAO);
  glBindVertexArray(VAO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  // glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
  // 1. then set the vertex attributes pointers
  glBufferData(GL_ARRAY_BUFFER, mesh->vertexArraySize * sizeof(float),
               mesh->vertexArray, GL_STATIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
  // glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float),
  //                       (void *)(3 * sizeof(float)));
  glEnableVertexAttribArray(0);
  // glEnableVertexAttribArray(1);
  // 2. use our shader program when we want to render an object
  glUseProgram(shaderProgram);

  int height;
  int width;
  glfwGetFramebufferSize(window, &width, &height);
  glViewport(0, 0, width, height);
  while (!glfwWindowShouldClose(window)) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    processInput(window);

    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glUseProgram(shaderProgram);
    glBindVertexArray(VAO);
    glDrawArrays(GL_TRIANGLES, 0, mesh->numTriangles*3);

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
