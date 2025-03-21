#include "../include/problem.h"

void fullSystemConstrain(fullSystem *mySystem, int myNode, double myValue) {
  double **A, *B;
  int i, size;

  A = mySystem->A;
  B = mySystem->B;
  size = mySystem->size;

  for (i = 0; i < size; i++) {
    B[i] -= myValue * A[i][myNode];
    A[i][myNode] = 0;
  }

  for (i = 0; i < size; i++)
    A[myNode][i] = 0;

  A[myNode][myNode] = 1;
  B[myNode] = myValue;
}


void elasticityAddBoundaryCondition(problem *theProblem, char *nameDomain, boundaryType type, double value) {
  int iDomain = geoGetDomain(theProblem->geometry, nameDomain);
  if (iDomain == -1) {
    Error("Undefined domain :-(");
  }

  boundaryCondition* theBoundary = malloc(sizeof(boundaryCondition));
  theBoundary->domain = theProblem->geometry->theDomains[iDomain];
  theBoundary->value = value;
  theBoundary->type = type;
  theProblem->nBoundaryConditions++;
  int size = theProblem->nBoundaryConditions;

  if (theProblem->conditions == NULL)
      theProblem->conditions = malloc(size*sizeof(boundaryCondition*));
  else
      theProblem->conditions = realloc(theProblem->conditions, size*sizeof(boundaryCondition*));
  theProblem->conditions[size-1] = theBoundary;

  int shift=-1;
  if (type == DIRICHLET_X)  shift = 0;
  if (type == DIRICHLET_Y)  shift = 1;
  if (shift == -1) return;
  int *elem = theBoundary->domain->elem;
  int nElem = theBoundary->domain->nElem;
  for (int e=0; e<nElem; e++) {
      for (int i=0; i<2; i++) {
          int node = theBoundary->domain->mesh->elem[2*elem[e]+i];
          theProblem->constrainedNodes[2*node+shift] = size-1; }}
}
