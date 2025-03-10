#pragma  once
#include <gmshc.h>
#include <stdlib.h>
typedef struct {
    float x;
    float y;
    float z;  // Adding z coordinate (set to 0.0 for 2D meshes)
} Vertex;

typedef struct {
    int v1;
    int v2;
    int v3;
} Triangle;

typedef struct {
    int numVertices;
    int numTriangles;
    float* vertexArray;  // Flattened array for OpenGL
    int vertexArraySize; // Size in number of floats
} Mesh;

Mesh* readMesh(const char* filename);
void freeMesh(Mesh* mesh);


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
} MeshSettings;

Mesh* generate_mesh(MeshSettings* s);
void normalize_mesh(Mesh* mesh);
