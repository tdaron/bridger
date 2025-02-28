#pragma once
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "../mesh/mesh.h"
GLFWwindow* createWindow(int width, int height, char* name);
GLuint compile_shader(const char *filename, GLenum shader_type);
unsigned int load_mesh_into_vao(unsigned int VAO, Mesh* mesh);
unsigned int get_program(char* vertex_shader_file, char* fragment_shader_file);
unsigned int create_vao();
