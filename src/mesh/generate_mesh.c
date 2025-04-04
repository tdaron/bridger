#include "gmshc.h"
#include "mesh.h"
#include <log.h>
#include <math.h>
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

// New function that only identifies tank edges
void find_tank_edges(Mesh *mesh, int **tankEdges, int *ntankEdges,
                     MeshSettings *s) {
  double maxY = -1e12;
  size_t nNode = mesh->nNode;
  size_t nEdge = mesh->nEdge;
  size_t *edgeNode = mesh->edgeNode;
  double *xyz = mesh->xyz;
  // Find max Y coordinate
  for (size_t i = 0; i < nNode; i++) {
    if (xyz[3 * i + 1] > maxY)
      maxY = xyz[3 * i + 1];
  }

  // Allocate memory for tank edges
  *tankEdges = malloc(sizeof(int) * nEdge);
  *ntankEdges = 0;

  // Identify tank edges
  for (int i = 0; i < nEdge; i++) {
    int n1 = edgeNode[2 * i] - 1;
    int n2 = edgeNode[2 * i + 1] - 1;
    double x1 = xyz[3 * n1];
    double y1 = xyz[3 * n1 + 1];
    double x2 = xyz[3 * n2];
    double y2 = xyz[3 * n2 + 1];

    // Identify tank edges
    if (fabs(y1 - maxY) <= 1e-5 && fabs(y2 - maxY) <= 1e-5) {
      if ((x1 > s->tankX && x1 < s->tankX + s->tankLength) ||
          (x2 > s->tankX && x2 < s->tankX + s->tankLength)) {
        (*tankEdges)[*ntankEdges] = i;
        (*ntankEdges)++;
      }
    }
  }
}

// Modified original function that no longer handles tank edges
void find_pillars_and_extremities(double *xyz, size_t nNode, size_t *edgeNode,
                                  size_t nEdge, int ***pillarEdges,
                                  int **nPillarEdges, int ***extremityEdges,
                                  int **nExtremityEdges, MeshSettings *s) {
  double minY = 1e12;
  double maxY = -1e12;
  double minX = 1e12;
  double maxX = -1e12;

  // Find min and max Y coordinates
  for (size_t i = 0; i < nNode; i++) {
    if (xyz[3 * i + 1] < minY)
      minY = xyz[3 * i + 1];
    if (xyz[3 * i + 1] > maxY)
      maxY = xyz[3 * i + 1];
  }

  // Find min and max X coordinates
  for (size_t i = 0; i < nNode; i++) {
    if (xyz[3 * i] < minX)
      minX = xyz[3 * i];
    if (xyz[3 * i] > maxX)
      maxX = xyz[3 * i];
  }

  // Identify pillar edges and extremity edges
  *pillarEdges = malloc(s->pillarsNumber * sizeof(int *));
  *nPillarEdges = malloc(s->pillarsNumber * sizeof(int));
  for (int i = 0; i < s->pillarsNumber; i++) {
    (*pillarEdges)[i] = malloc(nEdge * sizeof(int));
    (*nPillarEdges)[i] = 0;
  }

  *extremityEdges = malloc(2 * sizeof(int *));
  (*extremityEdges)[0] = malloc(nEdge * sizeof(int));
  (*extremityEdges)[1] = malloc(nEdge * sizeof(int));
  *nExtremityEdges = malloc(sizeof(int) * 2);
  (*nExtremityEdges)[0] = 0;
  (*nExtremityEdges)[1] = 0;

  for (int i = 0; i < nEdge; i++) {
    int n1 = edgeNode[2 * i] - 1;
    int n2 = edgeNode[2 * i + 1] - 1;
    double x1 = xyz[3 * n1];
    double y1 = xyz[3 * n1 + 1];
    double x2 = xyz[3 * n2];
    double y2 = xyz[3 * n2 + 1];

    // Identify pillar edges
    if (fabs(y1 - minY) < 1 && fabs(y2 - minY) < 1) {
      for (int j = 0; j < s->pillarsNumber; j++) {
        double px = get_pillar_x(j, s);
        if ((x1 >= px && x1 <= px + s->pillarsWidth) ||
            (x2 >= px && x2 <= px + s->pillarsWidth)) {
          (*pillarEdges)[j][(*nPillarEdges)[j]] = i;
          (*nPillarEdges)[j]++;
          break;
        }
      }
    }

    // Identify extremity edges
    if ((fabs(x1 - minX) < 1 && fabs(x2 - minX) < 1)) {
      if (fabs(y1 - maxY) <= 1e-5 && fabs(y2 - maxY) <= 1e-5) {
        continue;
      }
      (*extremityEdges)[0][(*nExtremityEdges)[0]] = i;
      (*nExtremityEdges)[0]++;
    }
    if ((fabs(x1 - maxX) < 1 && fabs(x2 - maxX) < 1)) {
      if (fabs(y1 - maxY) <= 1e-5 && fabs(y2 - maxY) <= 1e-5) {
        continue;
      }
      (*extremityEdges)[1][(*nExtremityEdges)[1]] = i;
      (*nExtremityEdges)[1]++;
    }
  }
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
    int idHole = gmshModelOccAddDisk(s->holeX, s->holeY, 0, s->holeRadius, s->holeRadius, -1, NULL,
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
    fprintf(file, "%6d : %14.7e %14.7e\n", i, xyz[3 * node[i] - 3],
            xyz[3 * node[i] - 2]);
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

  // Analyze geometry to create domains
  int **pillarEdges, **extremityEdges;
  int *nPillarEdges, *nExtremityEdges;
  size_t nIndNode = nNode;
  find_pillars_and_extremities(xyz, nNode, edgeNode, nEdge, &pillarEdges,
                               &nPillarEdges, &extremityEdges, &nExtremityEdges,
                               global_s);

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

  // Write domains to file
  int totalDomains =
      global_s->pillarsNumber + 2 + 1; // Number of pillars + 2 extremities
  fprintf(file, "Number of domains %d\n", totalDomains);

  // Write pillar domains
  for (int i = 0; i < global_s->pillarsNumber; i++) {
    fprintf(file, "  Domain : %6d\n", i + 1);
    fprintf(file, "  Name : Pillar%d\n", i + 1);
    fprintf(file, "  Number of elements : %6u\n", nPillarEdges[i]);
    for (int j = 0; j < nPillarEdges[i]; j++) {
      fprintf(file, "%6d ", pillarEdges[i][j]);
      if ((j + 1) % 10 == 0) {
        fprintf(file, "\n");
      }
    }
    fprintf(file, "\n");
    free(pillarEdges[i]);
  }

  // Write extremity domains
  for (int i = 0; i < 2; i++) {
    fprintf(file, "  Domain : %6d \n", global_s->pillarsNumber + i + 1);
    fprintf(file, "  Name : Extremity%d\n", i);
    fprintf(file, "  Number of elements : %6u\n", nExtremityEdges[i]);
    for (int j = 0; j < nExtremityEdges[i]; j++) {
      fprintf(file, "%6d ", extremityEdges[i][j]);
      if ((j + 1) % 10 == 0) {
        fprintf(file, "\n");
      }
    }
    fprintf(file, "\n");
    free(extremityEdges[i]);
  }

  free(pillarEdges);
  free(extremityEdges);
  free(nPillarEdges);

  // Create and return the mesh
  Mesh *mesh = malloc(sizeof(Mesh));
  mesh->numTriangles = nElem;
  mesh->numVertices = nNode;
  mesh->vertexArraySize =
      3 * 3 * nElem; // Each triangle has 3 vertices with 3 coordinates
  mesh->vertexArray = malloc(sizeof(float) * mesh->vertexArraySize);
  mesh->xyz = xyz;
  mesh->edgeNode = edgeNode;
  mesh->nEdge = nEdge;
  mesh->nNode = nIndNode;
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

  // Write tank domain

  int *tankEdges;
  int nTankEdges;
  find_tank_edges(mesh, &tankEdges, &nTankEdges, global_s);
  fprintf(file, "  Domain : %6d \n", global_s->pillarsNumber + 2 + 1);
  fprintf(file, "  Name : Tank\n");
  fprintf(file, "  Number of elements : %6u\n", nTankEdges);
  for (int j = 0; j < nTankEdges; j++) {
    fprintf(file, "%6d ", tankEdges[j]);
    if ((j + 1) % 10 == 0) {
      fprintf(file, "\n");
    }
  }
  fclose(file);


  gmshFree(edgeElem);
  // gmshFree(edgeNode);
  // gmshFree(xyz);
  gmshFree(trash);
  gmshFree(node);
  gmshFree(elem);

  *scale = normalize_mesh(mesh);
  return mesh;
}
