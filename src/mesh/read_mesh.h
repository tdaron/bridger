#pragma  once
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
    Vertex* vertices;
    int numVertices;
    Triangle* triangles;
    int numTriangles;
    float* vertexArray;  // Flattened array for OpenGL
    int vertexArraySize; // Size in number of floats
} Mesh;

Mesh* readMeshForOpenGL(const char* filename);
void freeMesh(Mesh* mesh);
