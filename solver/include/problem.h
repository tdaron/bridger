#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "shape.h"
#include "integration.h"

#define MAXNAME 256

typedef enum { FEM_TRIANGLE, FEM_QUAD, FEM_EDGE } elementType;
typedef enum { PLANAR_STRESS, PLANAR_STRAIN, AXISYM } elasticCase;
typedef enum { DIRICHLET_X, DIRICHLET_Y, NEUMANN_X, NEUMANN_Y } boundaryType;

#define ErrorScan(a) errorScan(a, __LINE__, __FILE__)
#define Error(a) error(a, __LINE__, __FILE__)

/**
 * @brief Structure containing the nodes of the global mesh
 * @param nNodes Number of nodes
 * @param X x-coordinates of the nodes
 * @param Y y-coordinates of the nodes
 */
typedef struct {
  int nNodes;
  double *X;
  double *Y;
} nodes;


/**
 * @brief Structure containing the elements of a mesh. Typically there are 2 meshes in a project:
 * one for the edges (segments) and one for the elements inside (e.g. triangles or quads)
 * @param nLocalNode Number of nodes per element, e.g. 2 for a segment (e.g. elements are edges), 3 for a triangle, 4 for a quad
 * @param nElem Number of elements
 * @param elem List of elements' nodes indices (size: `nElem * nLocalNode`).
 * It is arranged as follows: `elem[i*nLocalNode + j]` is the index of the j-th node of the i-th element.
 * For example, for a triangle, `elem[3*i]`, `elem[3*i+1]` and `elem[3*i+2]` are the indices of the 3 nodes of the i-th triangle
 * @param nodes Pointer to the nodes' structure
 */
typedef struct {
  int nLocalNode;
  int nElem;
  int *elem;
  nodes *nodes;
} mesh;

/**
 * @brief Structure containing one domain (a boundary) of the general geometry (e.g. the left side of a rectangle)
 * @param mesh Pointer to the mesh structure: contains elements of 2 nodes (edges)
 * @param nElem Number of elements in the domain
 * @param elem List of domain's elements' indices
 * @param name Name of the domain
 */
typedef struct {
  mesh *mesh;
  int nElem;
  int *elem;
  char name[MAXNAME];
} domain;

typedef struct {
  domain *domain;
  boundaryType type;
  double value;
} boundaryCondition;

/**
 * @brief Structure containing the general geometry of the problem
 * @param elementType Type of elements used in the mesh (e.g. triangles or quads)
 * @param theNodes Pointer to the nodes structure containing the nodes of the global mesh
 * @param theElements Pointer to the elements structure
 * @param theEdges Pointer to the edges structure
 * @param nDomains Number of domains in the geometry
 * @param theDomains List of domains
 */
typedef struct {
  elementType elementType;
  nodes *theNodes;
  mesh *theElements;
  mesh *theEdges;
  int nDomains;
  domain **theDomains;
} geo;

/**
 * @brief Structure containing the full system of equations
 * @param B Right-hand side of the system
 * @param A Stiffness matrix of the system
 * @param size Size of the system ( = 2*Nnodes)
 */
typedef struct {
  double *B;
  double **A;
  int size;
  int band;
} linearSystem;

typedef struct {
  double E, nu, rho, g;
  double A, B, C;
  int planarStrainStress;
  int nBoundaryConditions;
  boundaryCondition **conditions;
  int *constrainedNodes;
  geo *geometry;
  discrete *space;
  integration *rule;
  linearSystem *system;

  // New
  double* soluce;
  double* residuals;

  discrete *spaceEdge;
  integration *ruleEdge;

  int* renumOld2New;
  int* renumNew2Old;

} problem;

void geoInitialize(geo *theGeometry);
void geoFree(geo *theGeometry);
void geoFree(geo *theGeometry);
void geoSetDomainName(geo *theGeometry, int iDomain, char *name);
int geoGetDomain(geo *theGeometry, char *name);
void geoMeshPrint(geo *theGeometry);

linearSystem *femFullSystemCreate(int size);
linearSystem *bandFullSystemCreate(int size);
void fullSystemFree(linearSystem *theSystem);
void fullSystemAlloc(linearSystem *mySystem, int size);
void bandSystemAlloc(linearSystem *mySystem, int size);
void fullSystemInit(linearSystem *mySystem);
void bandSystemInit(linearSystem *mySystem);
void fullSystemPrint(linearSystem *mySystem);

problem *elasticityCreate(geo *theGeometry, double E, double nu, double rho, double g, elasticCase iCase, int makeBanded);
void elasticityFree(problem *theProblem);
void elasticityPrint(problem *theProblem);

void error(char *text, int line, char *file);
void errorScan(int test, int line, char *file);

double vecMin(double *x, int n);
double vecMax(double *x, int n);

void fullSystemConstrain(linearSystem *mySystem, int myNode, double myValue);
void bandSystemConstrain(linearSystem *mySystem, int myNode, double myValue);
void elasticityAddBoundaryCondition(problem *theProblem, char *nameDomain, boundaryType type, double value);

discrete *discreteCreate(int n, elementType type);

double elasticityIntegrate(problem *theProblem, double (*f)(double x, double y));
integration *integrationCreate(int n, elementType type);
void integrationFree(integration *theRule);

geo *geoMeshRead(const char *filename);
void solutionWrite(int nNodes, int nfields, double *data, const char *filename);
int solutiondRead(int allocated_size, double *value, const char *filename);

void elasticityAssembleElements(problem *theProblem);
void elasticityAssembleNeumann(problem *theProblem);
void elasticityAssembleElementsBand(problem *theProblem);

double *fullSystemEliminate(linearSystem *mySystem);
double *bandSystemEliminate(linearSystem *myBand);
double *bandSystemCholesky(linearSystem *myBand);

double * elasticityForces(problem *theProblem);
double *elasticitySolve(problem *theProblem, int makeBanded);

// Renumbering related functions
void problemXRenumber(problem *theProblem);
void rcm(problem* theProblem);
