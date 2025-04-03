#include "gmshc.h"
#include "mesh.h"
#include <log.h>
#include <stdlib.h>
Mesh *load_mesh_and_write_to_file(double *scale);

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

Mesh *generate_mesh(MeshSettings *s, double *scale) {
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
    int idHole = gmshModelOccAddDisk(s->holeX, s->holeY, 0, 0.6, 0.6, -1, NULL,
                                     0, NULL, 0, &ierr);
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

  return load_mesh_and_write_to_file(scale);
}

#define FILENAME "./solver/data/mesh.txt"

Mesh *load_mesh_and_write_to_file(double *scale) {
  FILE *file = fopen(FILENAME, "w");
  int ierr;

  // Getting vertices
  size_t nNode, n, m, *node;
  double *xyz, *trash;
  gmshModelMeshGetNodes(&node, &nNode, &xyz, &n, &trash, &m, -1, -1, 0, 0,
                        &ierr);
  ErrorGmsh(ierr);

  // Write nodes to file
  fprintf(file, "Number of nodes %zu \n", nNode);
  for (int i = 0; i < nNode; i++) {
    fprintf(file, "%6d : %14.7e %14.7e %14.7e\n", i, xyz[3 * node[i] - 3],
            xyz[3 * node[i] - 2], xyz[3 * node[i] - 1]);
  }

  // Getting edges
  size_t nEdge, *edgeElem;
  size_t *edgeNode;
  gmshModelMeshGetElementsByType(1, &edgeElem, &nEdge, &edgeNode, &n, -1, 0, 1,
                                 &ierr);
  ErrorGmsh(ierr);

  // Write edges to file
  fprintf(file, "Number of edges %zu \n", nEdge);
  for (int i = 0; i < nEdge; i++) {
    fprintf(file, "%6d : %6d %6d \n", i, (int)edgeNode[2 * i] - 1,
            (int)edgeNode[2 * i + 1] - 1);
  }
  gmshFree(edgeElem);
  gmshFree(edgeNode);

  // Getting triangles
  size_t nElem, *elem;
  gmshModelMeshGetElementsByType(2, &elem, &nElem, &node, &nNode, -1, 0, 1,
                                 &ierr);
  ErrorGmsh(ierr);

  // Write triangles to file
  fprintf(file, "Number of triangles %zu \n", nElem);
  for (int i = 0; i < nElem; i++) {
    fprintf(file, "%6d : %6d %6d %6d\n", i, (int)node[3 * i] - 1,
            (int)node[3 * i + 1] - 1, (int)node[3 * i + 2] - 1);
  }

  // Getting domains (1D entities)
  int *dimTags;
  size_t nDomains;
  gmshModelGetEntities(&dimTags, &nDomains, 1, &ierr);
  ErrorGmsh(ierr);
  nDomains = nDomains / 2;

  // Write domains to file
  fprintf(file, "Number of domains %zu\n", nDomains);

  for (int i = 0; i < nDomains; i++) {
    int dim = dimTags[2 * i + 0];
    int tag = dimTags[2 * i + 1];

    // Get domain elements
    int *elementType;
    size_t nElementType, **elementTags, *nElementTags, nnElementTags,
        **nodesTags, *nNodesTags, nnNodesTags;
    gmshModelMeshGetElements(&elementType, &nElementType, &elementTags,
                             &nElementTags, &nnElementTags, &nodesTags,
                             &nNodesTags, &nnNodesTags, dim, tag, &ierr);
    ErrorGmsh(ierr);

    fprintf(file, "  Domain : %6d \n", i);
    fprintf(file, "  Name : Entity %d \n", tag - 1);
    fprintf(file, "  Number of elements : %6zu\n", nElementTags[0]);

    // Write domain elements
    for (int j = 0; j < nElementTags[0]; j++) {
      fprintf(file, "%6zu",
              elementTags[0][j] - 1); // Adjust element indices to be 0-based
      if ((j + 1) != nElementTags[0] && (j + 1) % 10 == 0)
        fprintf(file, "\n");
    }
    fprintf(file, "\n");

    // Free allocated memory for domain elements
    gmshFree(nElementTags);
    gmshFree(nNodesTags);
    gmshFree(elementTags[0]);
    gmshFree(elementTags);
    gmshFree(nodesTags[0]);
    gmshFree(nodesTags);
    gmshFree(elementType);
  }
  gmshFree(dimTags);

  fclose(file);

  // Create and return the mesh
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

  *scale = normalize_mesh(mesh);
  return mesh;
}
