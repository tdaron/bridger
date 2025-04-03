#pragma once

/**
 * @brief Structure the parent element for a given element type (e.g. triangles or quads or edges)
 * @param n Number of nodes of an element (e.g. 3 for a triangle)
 * @param x2 Function that computes the parent coordinates of the element nodes (if 2D)
 * @param phi2 Function that computes the shape functions at a given integration point (if 2D)
 * @param dphi2dx Function that computes the derivatives of the shape functions at a given integration point (if 2D)
 * @param x Function that computes the parent coordinates of the element nodes (if 1D)
 * @param phi Function that computes the shape functions at a given integration point (if 1D)
 * @param dphidx Function that computes the derivatives of the shape functions at a given integration point (if 1D)
 */
typedef struct {
  int n;
  void (*x2)(double *xsi, double *eta);
  void (*phi2)(double xsi, double eta, double *phi);
  void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
  void (*x)(double *xsi);
  void (*phi)(double xsi, double *phi);
  void (*dphidx)(double xsi, double *dphidxsi);
} discrete;

void _q1c0_x(double *xsi, double *eta);
void _q1c0_phi(double xsi, double eta, double *phi);
void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta);

void _p1c0_x(double *xsi, double *eta);
void _p1c0_phi(double xsi, double eta, double *phi);
void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta);

void _e1c0_x(double *xsi);
void _e1c0_phi(double xsi, double *phi);
void _e1c0_dphidx(double xsi, double *dphidxsi);

void discreteFree(discrete *theSpace);

void discreteXsi(discrete *mySpace, double *xsi);
void discretePhi(discrete *mySpace, double xsi, double *phi);

void discreteXsi2(discrete *mySpace, double *xsi, double *eta);
void discretePhi2(discrete *mySpace, double xsi, double eta, double *phi);
void discreteDphi2(discrete *mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);

void discretePrint(discrete *mySpace);
