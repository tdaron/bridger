#pragma once
#include <gmshc.h>
typedef struct {
  float x;
  float y;
  float z; // Adding z coordinate (set to 0.0 for 2D meshes)
} Vertex;

typedef struct {
  int v1;
  int v2;
  int v3;
} Triangle;

typedef struct {
  int numVertices;
  int numTriangles;
  float *vertexArray;  // Flattened array for OpenGL
  int vertexArraySize; // Size in number of floats
  float centers[3];
  float scale;
  double *xyz;
  size_t nEdge;
  size_t *edgeNode;
  size_t nNode;
} Mesh;

Mesh *readMesh(const char *filename);
void freeMesh(Mesh *mesh);

typedef struct {
  int pillarsNumber;
  double pillarsWidth;
  double pillarsHeight;
  double bridgeWidth;
  double bridgeHeight;
  double offset;
  double baseElementSize;
  double preciseElementSize;
  double precisionRadius;
  double holeX;
  double holeY;
  double tankLength;
  double tankWeight;
  double tankX;
  double holeRadius;
} MeshSettings;

Mesh *generate_mesh(MeshSettings *s, double *scale);
double normalize_mesh(Mesh *mesh);
double getSize(double x, double y, MeshSettings *s);
double dn(double coord, Mesh *mesh, int k);
double *compute_field(Mesh *mesh, MeshSettings *s);
void find_tank_edges(Mesh *mesh, int **tankEdges, int *ntankEdges,
                     MeshSettings *s);
