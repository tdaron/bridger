#include "gmshc.h"
#include "mesh.h"
#include <log.h>
#include <stdlib.h>
Mesh *load_mesh();

void femErrorGmsh(int ierr, int line, char *file) {
  if (ierr == 0)
    return;
  log_error(
      "\n-------------------------------------------------------------------"
      "------------- ");
  log_error("\n  Error in %s at line %d : \n  error code returned by gmsh %d\n",
            file, line, ierr);
  log_error(
      "---------------------------------------------------------------------"
      " Yek Yek !! \n\n");
  gmshFinalize(NULL);
}

#define ErrorGmsh(a) femErrorGmsh(a, __LINE__, __FILE__)

double get_pillar_x(int i, MeshSettings *s) {
  return (double)(s->bridgeWidth - s->pillarsWidth - s->offset) * i /
             (s->pillarsNumber - 1) +
         s->offset / 2;
}

double getSize(double x, double y, MeshSettings *s) {
  double size = s->baseElementSize;

  double corners[4][2] = {{0, 0},
                          {0, s->pillarsHeight},
                          {s->pillarsWidth, 0},
                          {s->pillarsWidth, s->pillarsHeight}};

  for (int i = 0; i < s->pillarsNumber; i++) {
    double px = get_pillar_x(i, s);
    for (int j = 0; j < 4; j++) {
      // for each corner
      double cx = px + corners[j][0];
      double cy = -s->pillarsHeight + corners[j][1];

      double dx = x - cx;
      double dy = y - cy;

      double dist = dx * dx + dy * dy;
      if (dist < s->precisionRadius) {
        double new_size = dist * (s->baseElementSize - s->preciseElementSize) /
                              s->precisionRadius +
                          s->preciseElementSize;
        if (new_size < size)
          size = new_size;
      }
    }
  }

  if (x < s->precisionRadius) {
    double new_size =
        x * (s->baseElementSize - s->preciseElementSize) / s->precisionRadius +
        s->preciseElementSize;
    if (new_size < size)
      size = new_size;
  }
  if (x > s->bridgeWidth - s->precisionRadius) {
    double new_size = (s->bridgeWidth - x) *
                          (s->baseElementSize - s->preciseElementSize) /
                          s->precisionRadius +
                      s->preciseElementSize;
    if (new_size < size)
      size = new_size;
  }

  return size;
}

// cursed. who cares. (might wanna use the data attr from meshSize instead.)
MeshSettings *global_s;
double meshSize(int dim, int tag, double x, double y, double z, double lc,
                void *data) {

  return getSize(x, y, global_s);
}

Mesh *generate_mesh(MeshSettings *s) {
  // cursed. who cares.
  global_s = s;
  int ierr;

  // Creating the model
  gmshModelAdd("Bridge", &ierr);
  ErrorGmsh(ierr);
  gmshModelMeshSetSizeCallback(meshSize, NULL, &ierr);
  ErrorGmsh(ierr);

  // Adding the main bridge part
  int idBridge = gmshModelOccAddRectangle(0, 0, 0, s->bridgeWidth,
                                          s->bridgeHeight, -1, 0, &ierr);
  int bridge[] = {2, idBridge};

  for (int i = 0; i < s->pillarsNumber; i++) {
    double x = get_pillar_x(i, s);
    int idPillar =
        gmshModelOccAddRectangle(x, -s->pillarsHeight, 0, s->pillarsWidth,
                                 s->pillarsHeight, -1, 0, &ierr);
    int pillar[] = {2, idPillar};
    gmshModelOccFuse(bridge, 2, pillar, 2, NULL, NULL, NULL, NULL, NULL, -1, 1,
                     1, &ierr);
  }

  if (s->holeX != 0) {
    int idHole = gmshModelOccAddDisk(s->holeX, s->holeY, 0, 0.6, 0.6, -1, NULL, 0,
                                     NULL, 0, &ierr);
    int hole[] = {2, idHole};
    gmshModelOccCut(bridge, 2, hole, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1,
                    &ierr);
  }

  // Synchronizing CAD representation with gmsh internal models
  gmshModelOccSynchronize(&ierr);
  gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);

  // meshing
  gmshModelMeshGenerate(2, &ierr);
  ErrorGmsh(ierr);

  return load_mesh();
}

Mesh *load_mesh() {
  int ierr;

  // Getting vertices
  size_t nNode, n, m, *node;
  double *xyz, *trash;
  gmshModelMeshGetNodes(&node, &nNode, &xyz, &n, &trash, &m, -1, -1, 0, 0,
                        &ierr);
  ErrorGmsh(ierr);

  // getting triangles
  size_t nElem, *elem;
  gmshModelMeshGetElementsByType(2, &elem, &nElem, &node, &nNode, -1, 0, 1,
                                 &ierr);
  ErrorGmsh(ierr);

  Mesh *mesh = malloc(sizeof(Mesh));
  mesh->numTriangles = nElem;
  mesh->numVertices = nNode;
  mesh->vertexArraySize =
      3 * 3 * nElem; // Each triangle has 3 vertices with 3 coordinates
  mesh->vertexArray = malloc(sizeof(float) * mesh->vertexArraySize);
  for (size_t i = 0; i < nElem; i++) {
    // For each triangle
    for (size_t j = 0; j < 3; j++) {
      // 3 vertices
      int v = node[i * 3 + j] - 1;

      for (size_t k = 0; k < 3; k++) {
        // 3 dimensions
        double coord = xyz[v * 3 + k];
        mesh->vertexArray[i * 3 * 3 + j * 3 + k] = coord;
      }
    }
  }

  gmshFree(xyz);
  gmshFree(trash);
  gmshFree(node);
  gmshFree(elem);

  normalize_mesh(mesh);
  return mesh;
}
