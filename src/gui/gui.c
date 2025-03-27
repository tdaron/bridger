// clang-format off
#include <glad/glad.h>
#include <GLFW/glfw3.h>
// clang-format on
#include <stdio.h>
#include <stdlib.h>

// Callback function for window resize events
void framebuffer_size_callback(GLFWwindow *window, int width, int height) {
  // Preserve the original aspect ratio by letterboxing/pillarboxing
  float originalAspect =
      800.0f / 600.0f; // Original window size (assuming 800x600)
  float currentAspect = (float)width / (float)height;

  if (currentAspect > originalAspect) {
    // Window is wider than original aspect, use height as limiting factor
    int viewportWidth = (int)(height * originalAspect);
    int xOffset = (width - viewportWidth) / 2;
    glViewport(xOffset, 0, viewportWidth, height);
  } else {
    // Window is taller than original aspect, use width as limiting factor
    int viewportHeight = (int)(width / originalAspect);
    int yOffset = (height - viewportHeight) / 2;
    glViewport(0, yOffset, width, viewportHeight);
  }
}

GLFWwindow *createWindow(int width, int height, char *name) {
  glfwInit();
  glfwWindowHint(GLFW_RESIZABLE, 1);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_SCALE_TO_MONITOR, GL_TRUE);

  GLFWwindow *window = glfwCreateWindow(width, height, name, NULL, NULL);
  if (window == NULL) {
    const char *description;
    int code = glfwGetError(&description);
    printf("Failed to create GLFW window. Error code: %d, Description: %s\n",
           code, description);
    glfwTerminate();
    return NULL;
  }

  glfwMakeContextCurrent(window);
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
    printf("Failed to initialize GLAD\n");
    return NULL;
  }

  glViewport(0, 0, width, height);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  // Resize if the OS scales the window for highdpi
  glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
  glfwGetFramebufferSize(window, &width, &height);
  glViewport(0, 0, width, height);

  return window;
}

unsigned int create_vbo() {
  unsigned int vbo;
  glGenBuffers(1, &vbo);
  return vbo;
}
unsigned int create_vao() {
  unsigned int VAO;
  glGenVertexArrays(1, &VAO);
  return VAO;
}
