#include "../include/problem.h"

void geoInitialize(geo *theGeometry) {
  theGeometry->theNodes = NULL;
  theGeometry->theElements = NULL;
  theGeometry->theEdges = NULL;
  theGeometry->nDomains = 0;
  theGeometry->theDomains = NULL;
}

void geoFree(geo *theGeometry) {
  if (theGeometry->theNodes) {
    free(theGeometry->theNodes->X);
    free(theGeometry->theNodes->Y);
    free(theGeometry->theNodes);
  }
  if (theGeometry->theElements) {
    free(theGeometry->theElements->elem);
    free(theGeometry->theElements);
  }
  if (theGeometry->theEdges) {
    free(theGeometry->theEdges->elem);
    free(theGeometry->theEdges);
  }
  for (int i = 0; i < theGeometry->nDomains; i++) {
    free(theGeometry->theDomains[i]->elem);
    free(theGeometry->theDomains[i]);
  }
  free(theGeometry->theDomains);
}

// void geoFree(geo *theGeometry){
//   free(theGeometry->theElements->elem);
//   free(theGeometry->theElements);
//   free(theGeometry->theEdges->elem);
//   free(theGeometry->theEdges);
//   free(theGeometry->theNodes->X);
//   free(theGeometry->theNodes->Y);
//   free(theGeometry->theNodes);
//   for(int i = 0; i < theGeometry->nDomains; i++){
//     domain *theDomain = theGeometry->theDomains[i];
//     free(theDomain->elem);
//     free(theDomain);
//   }
//   free(theGeometry->theDomains);
//   free(theGeometry);
// }

void geoSetDomainName(geo *theGeometry, int iDomain, char *name) {
  if (iDomain >= theGeometry->nDomains)
    Error("Illegal domain number");
  if (geoGetDomain(theGeometry, name) != -1)
    Error("Cannot use the same name for two domains");
  sprintf(theGeometry->theDomains[iDomain]->name, "%s", name);
}

int geoGetDomain(geo *theGeometry, char *name) {
  int theIndex = -1;
  int nDomains = theGeometry->nDomains;
  for (int iDomain = 0; iDomain < nDomains; iDomain++) {
    domain *theDomain = theGeometry->theDomains[iDomain];
    if (strncasecmp(name, theDomain->name, MAXNAME) == 0)
      theIndex = iDomain;
  }
  return theIndex;
}


void geoMeshPrint(geo *theGeometry) {
  nodes *theNodes = theGeometry->theNodes;
  if (theNodes != NULL) {
    printf("Number of nodes %d \n", theNodes->nNodes);
    for (int i = 0; i < theNodes->nNodes; i++) {
      printf("%6d : %14.7e %14.7e \n", i, theNodes->X[i], theNodes->Y[i]);
    }
  }
  mesh *theEdges = theGeometry->theEdges;
  if (theEdges != NULL) {
    printf("Number of edges %d \n", theEdges->nElem);
    int *elem = theEdges->elem;
    for (int i = 0; i < theEdges->nElem; i++) {
      printf("%6d : %6d %6d \n", i, elem[2 * i], elem[2 * i + 1]);
    }
  }
  mesh *theElements = theGeometry->theElements;
  if (theElements != NULL) {
    if (theElements->nLocalNode == 3) {
      printf("Number of triangles %d \n", theElements->nElem);
      int *elem = theElements->elem;
      for (int i = 0; i < theElements->nElem; i++) {
        printf("%6d : %6d %6d %6d\n", i, elem[3 * i], elem[3 * i + 1], elem[3 * i + 2]);
      }
    }
    if (theElements->nLocalNode == 4) {
      printf("Number of quads %d \n", theElements->nElem);
      int *elem = theElements->elem;
      for (int i = 0; i < theElements->nElem; i++) {
        printf("%6d : %6d %6d %6d %6d\n", i, elem[4 * i], elem[4 * i + 1], elem[4 * i + 2], elem[4 * i + 3]);
      }
    }
  }
  int nDomains = theGeometry->nDomains;
  printf("Number of domains %d\n", nDomains);
  for (int iDomain = 0; iDomain < nDomains; iDomain++) {
    domain *theDomain = theGeometry->theDomains[iDomain];
    printf("  Domain : %6d \n", iDomain);
    printf("  Name : %s\n", theDomain->name);
    printf("  Number of elements : %6d\n", theDomain->nElem);
    for (int i = 0; i < theDomain->nElem; i++) {
      //         if (i != theDomain->nElem  && (i % 10) != 0)  printf(" - ");
      printf("%6d", theDomain->elem[i]);
      if ((i + 1) != theDomain->nElem && (i + 1) % 10 == 0)
        printf("\n");
    }
    printf("\n");
  }
}


linearSystem *fullSystemCreate(int size) {
  linearSystem *theSystem = malloc(sizeof(linearSystem));
  theSystem->band = -1;
  fullSystemAlloc(theSystem, size);
  fullSystemInit(theSystem);

  return theSystem;
}

linearSystem *bandSystemCreate(int size, int band) {
  linearSystem *theSystem = malloc(sizeof(linearSystem));
  theSystem->band = band;
  bandSystemAlloc(theSystem, size);
  bandSystemInit(theSystem);

  return theSystem;
}

void fullSystemFree(linearSystem *theSystem) {
  free(theSystem->A);
  free(theSystem->B);
  free(theSystem);
}

void fullSystemAlloc(linearSystem *mySystem, int size) {
  int i;
  double *elem = malloc(sizeof(double) * size * (size + 1));
  mySystem->A = malloc(sizeof(double *) * size);
  mySystem->B = elem;
  mySystem->A[0] = elem + size;
  mySystem->size = size;
  for (i = 1; i < size; i++)
    mySystem->A[i] = mySystem->A[i - 1] + size;
}

void bandSystemAlloc(linearSystem *mySystem, int size) {
  int i;
  int band = mySystem->band;
  if (band == -1) {
    fprintf(stderr, "Bandwidth set to -1 in bandSystemAlloc? What the fuck are you doing?\n");
    exit(1);
  }
  double *elem = malloc(sizeof(double) * size * (band + 1 + 1));
  mySystem->A = malloc(sizeof(double *) * size);
  mySystem->B = elem;
  mySystem->A[0] = elem + size;
  mySystem->size = size;
  for (i = 1; i < size; i++)
    mySystem->A[i] = mySystem->A[i - 1] + band;
}

void fullSystemInit(linearSystem *mySystem) {
  int i, size = mySystem->size;
  for (i = 0; i < size * (size + 1); i++)
    mySystem->B[i] = 0;
}

void bandSystemInit(linearSystem *mySystem) {
  int i;
  int size = mySystem->size;
  for (i = 0; i < size * (mySystem->band + 1 + 1); i++)
    mySystem->B[i] = 0.0;
}

void fullSystemPrint(linearSystem *mySystem) {
  double **A, *B;
  int i, j, size;

  A = mySystem->A;
  B = mySystem->B;
  size = mySystem->size;

  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++)
      if (A[i][j] == 0)
        printf("         ");
      else
        printf(" %+.1e", A[i][j]);
    printf(" :  %+.1e \n", B[i]);
  }
}

problem *elasticityCreate(
    geo *theGeometry,
    double E, double nu, double rho, double g,
    elasticCase iCase,
    int makeBanded)
{
  problem *theProblem = malloc(sizeof(problem));
  theProblem->E = E;
  theProblem->nu = nu;
  theProblem->g = g;
  theProblem->rho = rho;

  int nLocalNode = theGeometry->theElements->nLocalNode;
  int nNodes = theGeometry->theNodes->nNodes;

  if (iCase == PLANAR_STRESS) {
    theProblem->A = E / (1 - nu * nu);
    theProblem->B = E * nu / (1 - nu * nu);
    theProblem->C = E / (2 * (1 + nu));
  } else if (iCase == PLANAR_STRAIN || iCase == AXISYM) {
    theProblem->A = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
    theProblem->B = E * nu / ((1 + nu) * (1 - 2 * nu));
    theProblem->C = E / (2 * (1 + nu));
  }

  theProblem->planarStrainStress = iCase;
  theProblem->nBoundaryConditions = 0;
  theProblem->conditions = NULL;

  int size = 2 * nNodes;
  theProblem->constrainedNodes = malloc(size * sizeof(int));
  theProblem->soluce = malloc(size * sizeof(double));
  theProblem->residuals = malloc(size * sizeof(double));
  for (int i = 0; i < size; i++) {
    theProblem->constrainedNodes[i] = -1;
    theProblem->soluce[i] = 0;
    theProblem->residuals[i] = 0;
  }

  theProblem->geometry = theGeometry;
  if (theGeometry->theElements->nLocalNode == 3) {
    theProblem->space = discreteCreate(3, FEM_TRIANGLE);
    theProblem->rule = integrationCreate(3, FEM_TRIANGLE);
  }
  if (theGeometry->theElements->nLocalNode == 4) {
    theProblem->space = discreteCreate(4, FEM_QUAD);
    theProblem->rule = integrationCreate(4, FEM_QUAD);
  }
  // New
  theProblem->spaceEdge    = discreteCreate(2,FEM_EDGE);
  theProblem->ruleEdge     = integrationCreate(2,FEM_EDGE);

  theProblem->renumOld2New = malloc(sizeof(int)*nNodes);
  theProblem->renumNew2Old = malloc(sizeof(int)*nNodes);

  // Node Renumbering
  problemXRenumber(theProblem);

  // Calculate bandwith
  int max = 0;
  for (int i = 0; i < theGeometry->theElements->nElem; i++) {
    for (int j = 0; j < nLocalNode; j++) {

      for (int k = 0; k < nLocalNode; k++) {
        int diff =
            abs(theProblem->renumOld2New[theGeometry->theElements->elem[i * nLocalNode + j]] -
            theProblem->renumOld2New[theGeometry->theElements->elem[i * nLocalNode + k]]);
        if (diff > max) {
          max = diff;
        }
      }
    }
  }


  // Need to double the element id distance because of the 2 dofs per node
  // and add 1 for the diagonal which goes from a single element when there
  // is a single degree of freedom per node to a 2x2 block when there are 2.
  // For example, the matrix will always be at least of bandwidth 1.
  // We'll also include the diagonal in the bandwidth.
  int bandwidth = 2 * max + 1;

  if (makeBanded)
    theProblem->system = bandSystemCreate(size, bandwidth);
  else
    theProblem->system = fullSystemCreate(size);

  return theProblem;
}

void elasticityFree(problem *theProblem) {
  fullSystemFree(theProblem->system);
  integrationFree(theProblem->rule);
  discreteFree(theProblem->space);
  for (int i = 0; i < theProblem->nBoundaryConditions; i++) {
    free(theProblem->conditions[i]);
  }
  free(theProblem->conditions);
  free(theProblem->constrainedNodes);
  free(theProblem);
}


void elasticityPrint(problem *theProblem)
{
    printf("\n\n ======================================================================================= \n\n");
    printf(" Linear elasticity problem \n");
    printf("   Young modulus   E   = %14.7e [N/m2]\n",theProblem->E);
    printf("   Poisson's ratio nu  = %14.7e [-]\n",theProblem->nu);
    printf("   Density         rho = %14.7e [kg/m3]\n",theProblem->rho);
    printf("   Gravity         g   = %14.7e [m/s2]\n",theProblem->g);

    if (theProblem->planarStrainStress == PLANAR_STRAIN)  printf("   Planar strains formulation \n");
    if (theProblem->planarStrainStress == PLANAR_STRESS)  printf("   Planar stresses formulation \n");
    if (theProblem->planarStrainStress == AXISYM)         printf("   Axisymmetric formulation \n");

    printf("   Boundary conditions : \n");
    for(int i=0; i < theProblem->nBoundaryConditions; i++) {
          boundaryCondition *theCondition = theProblem->conditions[i];
          double value = theCondition->value;
          printf("  %20s :",theCondition->domain->name);
          if (theCondition->type==DIRICHLET_X)  printf(" imposing %9.2e as the horizontal displacement  \n",value);
          if (theCondition->type==DIRICHLET_Y)  printf(" imposing %9.2e as the vertical displacement  \n",value);
          if (theCondition->type==NEUMANN_X)    printf(" imposing %9.2e as the horizontal force desnity \n",value);
          if (theCondition->type==NEUMANN_Y)    printf(" imposing %9.2e as the vertical force density \n",value);}
    printf(" ======================================================================================= \n\n");
}

void error(char *text, int line, char *file) {
  printf("\n-------------------------------------------------------------------------------- ");
  printf("\n  Error in %s:%d at line %d : \n  %s\n", file, line, line, text);
  printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
  exit(69);
}

void errorScan(int test, int line, char *file) {
  if (test >= 0)
    return;

  printf("\n-------------------------------------------------------------------------------- ");
  printf("\n  Error in fscanf or fgets in %s:%d at line %d : \n", file, line, line);
  printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
  exit(69);
}

double vecMin(double *x, int n) {
  double myMin = x[0];
  int i;
  for (i = 1; i < n; i++)
    myMin = fmin(myMin, x[i]);
  return myMin;
}

double vecMax(double *x, int n) {
  double myMax = x[0];
  int i;
  for (i = 1; i < n; i++)
    myMax = fmax(myMax, x[i]);
  return myMax;
}
