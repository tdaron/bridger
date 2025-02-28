#include "gui.h"
#include "stdio.h"
#include "stdlib.h"

/**
@returns The program created.
**/
unsigned int get_program(char* vertex_shader_file, char* fragment_shader_file) {
  unsigned int vertexShader =
      compile_shader(vertex_shader_file, GL_VERTEX_SHADER);
  unsigned int fragmentShader =
      compile_shader(fragment_shader_file, GL_FRAGMENT_SHADER);
  unsigned int shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, vertexShader);
  glAttachShader(shaderProgram, fragmentShader);
  glLinkProgram(shaderProgram);
  glDeleteShader(vertexShader);
  glDeleteShader(fragmentShader);
  return shaderProgram;
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
