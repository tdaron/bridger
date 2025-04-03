#pragma once
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "../mesh/mesh.h"
#include "problem.h"
GLFWwindow* createWindow(int width, int height, char* name);
GLuint compile_shader(const char *filename, GLenum shader_type);
unsigned int load_mesh_into_vao(unsigned int VAO, Mesh* mesh, unsigned int VBO);
unsigned int get_program(char* vertex_shader_file, char* fragment_shader_file);
unsigned int create_vao();
unsigned int create_vbo();
unsigned int load_field_into_vao(unsigned int VAO, double *field, int size, unsigned int VBO);
unsigned int load_soluce_into_vao(unsigned int VAO, double *soluce, geo* geometry, unsigned int VBO, double scale);

void DrawImage(char* path, float x, float y, float scale);
