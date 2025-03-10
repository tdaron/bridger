#include "mesh.h"
#include <math.h>
#include <stdlib.h>

void getColor(double value, int numberOfColors, float *r, float *g, float *b) {
  if (value > 1)
    value = 1;
  value = value * (numberOfColors);
  value = (int)(value - 0.00000001);
  value = value / (numberOfColors - 1);

  value = 1 - value;
  if (value < 0)
    value = 0;
  if (value > 1)
    value = 1;
  *r = 3.5 * (1 - value) * (1 - value);
  *g = (1 - value) * (value) * 3.5;
  *b = value * value;
}

double glScale(double minimum, double maximum, double value) {
  if (value < minimum)
    return 0;
  if (minimum == maximum)
    return minimum;
  return (value - minimum) / fabs(maximum - minimum);
}


float* compute_field(Mesh* mesh, MeshSettings* s) {
  float *field = malloc(mesh->numTriangles * 9 * sizeof(float));

  for (int i = 0; i < mesh->numTriangles; i++) {

#define size(i, j)                                                             \
  getSize(dn(mesh->vertexArray[i * 9 + j * 3], mesh, 0),                       \
          dn(mesh->vertexArray[i * 9 + j * 3 + 1], mesh, 1), s)
    for (int j = 0; j < 3; j++) {
      float r;
      float g;
      float b;
      getColor(glScale(s->preciseElementSize, s->baseElementSize,
                       size(i, j)),
               50, &r, &g, &b);
      field[i * 9 + j * 3] = r;
      field[i * 9 + j * 3 + 1] = g;
      field[i * 9 + j * 3 + 2] = b;
    }
  }
  return field;
}
