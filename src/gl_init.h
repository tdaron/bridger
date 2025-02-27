#pragma once
#include <GLFW/glfw3.h>
GLFWwindow* createWindow(int width, int height, char* name);
GLuint compile_shader(const char *filename, GLenum shader_type);
