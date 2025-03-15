
/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#ifndef _FEM_H_
#define _FEM_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ErrorScan(a) femErrorScan(a, __LINE__, __FILE__)
#define Error(a) femError(a, __LINE__, __FILE__)
#define Warning(a) femWarning(a, __LINE__, __FILE__)
#define FALSE 0
#define TRUE 1
#define MAXNAME 256

typedef enum { FEM_TRIANGLE, FEM_QUAD, FEM_EDGE} femElementType;
typedef enum { DIRICHLET_X, DIRICHLET_Y } femBoundaryType;
typedef enum { PLANAR_STRESS, PLANAR_STRAIN, AXISYM } femElasticCase;

typedef struct {
  int nNodes;
  double *X;
  double *Y;
} femNodes;

typedef struct {
  int nLocalNode;
  int nElem;
  int *elem;
  femNodes *nodes;
} femMesh;

typedef struct {
  femMesh *mesh;
  int nElem;
  int *elem;
  char name[MAXNAME];
} femDomain;

typedef struct {
  double LxPlate, LyPlate;
  double h;
  femElementType elementType;
  double (*geoSize)(double x, double y);
  femNodes *theNodes;
  femMesh *theElements;
  femMesh *theEdges;
  int nDomains;
  femDomain **theDomains;
} femGeo;

typedef struct {
  int n;
  void (*x2)(double *xsi, double *eta);
  void (*phi2)(double xsi, double eta, double *phi);
  void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
  void (*x)(double *xsi);
  void (*phi)(double xsi, double *phi);
  void (*dphidx)(double xsi, double *dphidxsi);
} femDiscrete;

typedef struct {
  int n;
  const double *xsi;
  const double *eta;
  const double *weight;
} femIntegration;

typedef struct {
  double *B;
  double **A;
  int size;
} femFullSystem;

typedef struct {
  femDomain *domain;
  femBoundaryType type;
  double value;
} femBoundaryCondition;

typedef struct {
  double E, nu, rho, g;
  double A, B, C;
  int planarStrainStress;
  int nBoundaryConditions;
  femBoundaryCondition **conditions;
  int *constrainedNodes;
  femGeo *geometry;
  femDiscrete *space;
  femIntegration *rule;
  femFullSystem *system;
} femProblem;

void geoInitialize(femGeo *theGeometry);
void geoMeshPrint(femGeo *theGeometry);
femGeo *geoMeshRead(const char *filename);
void geoSetDomainName(femGeo *theGeometry, int iDomain, char *name);
int geoGetDomain(femGeo *theGeometry, char *name);
void geoFinalize();

femProblem *femElasticityCreate(femGeo *theGeometry, double E, double nu, double rho, double g, femElasticCase iCase);
void femElasticityFree(femProblem *theProblem);
void femElasticityPrint(femProblem *theProblem);
void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value);
double *femElasticitySolve(femProblem *theProblem);

femIntegration *femIntegrationCreate(int n, femElementType type);
void femIntegrationFree(femIntegration *theRule);

femDiscrete *femDiscreteCreate(int n, femElementType type);
void femDiscreteFree(femDiscrete *mySpace);
void femDiscretePrint(femDiscrete *mySpace);
void femDiscreteXsi2(femDiscrete *mySpace, double *xsi, double *eta);
void femDiscretePhi2(femDiscrete *mySpace, double xsi, double eta, double *phi);
void femDiscreteDphi2(femDiscrete *mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);

femFullSystem *femFullSystemCreate(int size);
void femFullSystemFree(femFullSystem *mySystem);
void femFullSystemPrint(femFullSystem *mySystem);
void femFullSystemInit(femFullSystem *mySystem);
void femFullSystemAlloc(femFullSystem *mySystem, int size);
double *femFullSystemEliminate(femFullSystem *mySystem);
void femFullSystemConstrain(femFullSystem *mySystem, int myNode, double value);

double femMin(double *x, int n);
double femMax(double *x, int n);
void femError(char *text, int line, char *file);
void femErrorScan(int test, int line, char *file);
void femWarning(char *text, int line, char *file);

void femSolutionWrite(int nNodes, int nfields, double *data, const char *filename);

#endif
