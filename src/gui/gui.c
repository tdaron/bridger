// clang-format off
#include <glad/glad.h>
#include <GLFW/glfw3.h>
// clang-format on
#include <stdio.h>
#include <stdlib.h>


// Callback function for window resize events
void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
  // Preserve the original aspect ratio by letterboxing/pillarboxing
  float originalAspect = 800.0f / 600.0f;  // Original window size (assuming 800x600)
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

  glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
  return window;
}

// Function to read the file into a string
char *read_file(const char *filename) {
  FILE *file = fopen(filename, "rb");
  if (!file) {
    perror("Unable to open file");
    return NULL;
  }

  fseek(file, 0, SEEK_END);
  long length = ftell(file);
  fseek(file, 0, SEEK_SET);

  char *content = (char *)malloc(length + 1);
  if (!content) {
    perror("Unable to allocate memory for file content");
    fclose(file);
    return NULL;
  }

  fread(content, 1, length, file);
  content[length] = '\0'; // Null-terminate the string

  fclose(file);
  return content;
}

// Function to compile the shader from a file
GLuint compile_shader(const char *filename, GLenum shader_type) {
  // Read the shader source from the file
  char *shader_source = read_file(filename);
  if (!shader_source) {
    return 0; // Return 0 on failure
  }

  // Create the shader object
  GLuint shader = glCreateShader(shader_type);
  if (shader == 0) {
    fprintf(stderr, "Error creating shader\n");
    free(shader_source);
    return 0;
  }

  // Set the shader source code
  glShaderSource(shader, 1, (const char **)&shader_source, NULL);

  // Compile the shader
  glCompileShader(shader);

  // Check for compilation errors
  GLint success;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
  if (!success) {
    GLint log_length;
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &log_length);
    char *log = (char *)malloc(log_length);
    glGetShaderInfoLog(shader, log_length, &log_length, log);
    fprintf(stderr, "Shader compilation failed: %s\n", log);
    free(log);
    glDeleteShader(shader);
    free(shader_source);
    return 0; // Return 0 on failure
  }

  // Free the shader source after compiling
  free(shader_source);

  return shader; // Return the compiled shader
}

