#include "../include/problem.h"


geo *geoMeshRead(const char *filename) {
  geo *theGeometry = malloc(sizeof(geo));
  geoInitialize(theGeometry);

  FILE *file = fopen(filename, "r");
  if (!file) {
    printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename);
    exit(-1);
  }

  int trash, *elem;

  nodes *theNodes = malloc(sizeof(nodes));
  theGeometry->theNodes = theNodes;
  ErrorScan(fscanf(file, "Number of nodes %d \n", &theNodes->nNodes));
  theNodes->X = malloc(sizeof(double) * (theNodes->nNodes));
  theNodes->Y = malloc(sizeof(double) * (theNodes->nNodes));
  for (int i = 0; i < theNodes->nNodes; i++) {
    ErrorScan(fscanf(file, "%d : %le %le \n", &trash, &theNodes->X[i], &theNodes->Y[i]));
  }

  mesh *theEdges = malloc(sizeof(mesh));
  theGeometry->theEdges = theEdges;
  theEdges->nLocalNode = 2;
  theEdges->nodes = theNodes;
  ErrorScan(fscanf(file, "Number of edges %d \n", &theEdges->nElem));
  theEdges->elem = malloc(sizeof(int) * theEdges->nLocalNode * theEdges->nElem);
  for (int i = 0; i < theEdges->nElem; ++i) {
    elem = theEdges->elem;
    ErrorScan(fscanf(file, "%6d : %6d %6d \n", &trash, &elem[2 * i], &elem[2 * i + 1]));
  }

  mesh *theElements = malloc(sizeof(mesh));
  theGeometry->theElements = theElements;
  theElements->nLocalNode = 0;
  theElements->nodes = theNodes;
  char elementType[MAXNAME];
  ErrorScan(fscanf(file, "Number of %s %d \n", elementType, &theElements->nElem));
  if (strncasecmp(elementType, "triangles", MAXNAME) == 0) {
    theElements->nLocalNode = 3;
    theElements->elem = malloc(sizeof(int) * theElements->nLocalNode * theElements->nElem);
    for (int i = 0; i < theElements->nElem; ++i) {
      elem = theElements->elem;
      ErrorScan(fscanf(file, "%6d : %6d %6d %6d \n", &trash, &elem[3 * i], &elem[3 * i + 1], &elem[3 * i + 2]));
    }
  }
  if (strncasecmp(elementType, "quads", MAXNAME) == 0) {
    theElements->nLocalNode = 4;
    theElements->elem = malloc(sizeof(int) * theElements->nLocalNode * theElements->nElem);
    for (int i = 0; i < theElements->nElem; ++i) {
      elem = theElements->elem;
      ErrorScan(fscanf(file, "%6d : %6d %6d %6d %6d \n", &trash, &elem[4 * i], &elem[4 * i + 1], &elem[4 * i + 2], &elem[4 * i + 3]));
    }
  }

  ErrorScan(fscanf(file, "Number of domains %d\n", &theGeometry->nDomains));
  int nDomains = theGeometry->nDomains;
  theGeometry->theDomains = malloc(sizeof(domain *) * nDomains);
  for (int iDomain = 0; iDomain < nDomains; iDomain++) {
    domain *theDomain = malloc(sizeof(domain));
    theGeometry->theDomains[iDomain] = theDomain;
    theDomain->mesh = theEdges;
    ErrorScan(fscanf(file, "  Domain : %6d \n", &trash));
    ErrorScan(fscanf(file, "  Name : %[^\n]s \n", (char *)&theDomain->name));
    ErrorScan(fscanf(file, "  Number of elements : %6d\n", &theDomain->nElem));
    theDomain->elem = malloc(sizeof(int) * 2 * theDomain->nElem);
    for (int i = 0; i < theDomain->nElem; i++) {
      ErrorScan(fscanf(file, "%6d", &theDomain->elem[i]));
      if ((i + 1) != theDomain->nElem && (i + 1) % 10 == 0)
        ErrorScan(fscanf(file, "\n"));
    }
  }

  fclose(file);
  return theGeometry;
}


void solutionWrite(int nNodes, int nfields, double *data, const char *filename) {
  FILE *file = fopen(filename, "w");
  if (!file) {
    printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename);
    exit(-1);
  }
  fprintf(file, "Size %d,%d\n", nNodes, nfields);
  for (int i = 0; i < nNodes; i++) {
    for (int j = 0; j < nfields - 1; j++) {
      fprintf(file, "%.18le,", data[i * nfields + j]);
    }
    fprintf(file, "%.18le", data[i * nfields + nfields - 1]);
    fprintf(file, "\n");
  }
  fclose(file);
}

int solutiondRead(int allocated_size, double *value, const char *filename) {
  FILE *file = fopen(filename, "r");
  if (!file) {
    printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename);
    exit(-1);
  }
  int nNodes, nFields;
  ErrorScan(fscanf(file, "Size %d,%d\n", &nNodes, &nFields));
  if (nNodes * nFields > allocated_size) {
    printf("Error: allocated size is %d, but the solution file has %d nodes and %d fields", allocated_size, nNodes, nFields);
    Error("The allocated size is too small for solutiondRead");
  }
  for (int i = 0; i < nNodes; i++) {
    for (int j = 0; j < nFields; j++)
      ErrorScan(fscanf(file, "%le,", &value[i * nFields + j]));
    ErrorScan(fscanf(file, "\n"));
  }
  printf("Reading solution of shape (%d,%d)\n", nNodes, nFields);
  fclose(file);
  return nNodes * nFields;
}
