#define STB_IMAGE_IMPLEMENTATION
#include "../external/stb_image.h"
#include "gui.h"

typedef struct {
  char *path;
  int width;
  int height;
  unsigned int texture;
} Image;
extern int windowWidth, windowHeight;

int numberOfImages;
Image loaded[256];
unsigned int imgVbo = -1;
unsigned int imgVao = -1;
unsigned int shaderProgram = -1;
void DrawImage(char *path, float x, float y, float scale) {
  if (imgVao == -1) {
    imgVao = create_vao();
    imgVbo = create_vbo();
    glBindVertexArray(imgVao);
    glBindBuffer(GL_ARRAY_BUFFER, imgVbo);

    shaderProgram = get_program("src/gui/shaders/vertex_text.vert",
                                "src/gui/shaders/fragment_text.frag");
    glUseProgram(shaderProgram);
    glUniform1i(glGetUniformLocation(shaderProgram, "texture1"),
                0); // set it manually
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float),
                          (void *)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float),
                          (void *)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
  }
  for (int i = 0; i < numberOfImages; i++) {
    Image *img = loaded + i;

    if (img->path == path) {

      glBindVertexArray(imgVao);
      glBindBuffer(GL_ARRAY_BUFFER, imgVbo);

      float dx = ((float)img->width / windowWidth) * scale;
      float dy = ((float)img->height / windowHeight) * scale;

      float vertices[] = {
          // positions        / texture coords
          x + dx, y + dy, 0.0f, 1.0f, 1.0f, // top right
          x + dx, y,      0.0f, 1.0f, 0.0f, // bottom right
          x,      y,      0.0f, 0.0f, 0.0f, // bottom left
          x,      y + dy, 0.0f, 0.0f, 1.0f, // top left
          x + dx, y + dy, 0.0f, 1.0f, 1.0f, // top right
          x,      y,      0.0f, 0.0f, 0.0f, // bottom left
      };
      glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
      glUseProgram(shaderProgram);

      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, img->texture);

      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      glDrawArrays(GL_TRIANGLES, 0, 6);
      // draw already loaded image
      return;
    }
  }
  int width, height, nrChannels;
  stbi_set_flip_vertically_on_load(1);
  unsigned char *data = stbi_load(path, &width, &height, &nrChannels, 0);
  printf("Loaded image. Width: %d Height: %d\n", width, height);

  unsigned int texture;
  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexParameteri(
      GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
      GL_REPEAT); // set texture wrapping to GL_REPEAT (default wrapping method)
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  // set texture filtering parameters
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA,
               GL_UNSIGNED_BYTE, data);
  glGenerateMipmap(GL_TEXTURE_2D);
  stbi_image_free(data);

  Image i = {
      .path = path, .width = width, .height = height, .texture = texture};

  loaded[numberOfImages] = i;
  numberOfImages++;
}
