#include "../include/shape.h"
#include "../include/problem.h"

void _q1c0_x(double *xsi, double *eta) {
  xsi[0] = 1.0;
  eta[0] = 1.0;
  xsi[1] = -1.0;
  eta[1] = 1.0;
  xsi[2] = -1.0;
  eta[2] = -1.0;
  xsi[3] = 1.0;
  eta[3] = -1.0;
}

void _q1c0_phi(double xsi, double eta, double *phi) {
  phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;
  phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
  phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
  phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
}

void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta) {
  dphidxsi[0] = (1.0 + eta) / 4.0;
  dphidxsi[1] = -(1.0 + eta) / 4.0;
  dphidxsi[2] = -(1.0 - eta) / 4.0;
  dphidxsi[3] = (1.0 - eta) / 4.0;
  dphideta[0] = (1.0 + xsi) / 4.0;
  dphideta[1] = (1.0 - xsi) / 4.0;
  dphideta[2] = -(1.0 - xsi) / 4.0;
  dphideta[3] = -(1.0 + xsi) / 4.0;
}

void _p1c0_x(double *xsi, double *eta) {
  xsi[0] = 0.0;
  eta[0] = 0.0;
  xsi[1] = 1.0;
  eta[1] = 0.0;
  xsi[2] = 0.0;
  eta[2] = 1.0;
}

void _p1c0_phi(double xsi, double eta, double *phi) {
  phi[0] = 1 - xsi - eta;
  phi[1] = xsi;
  phi[2] = eta;
}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta) {
  dphidxsi[0] = -1.0;
  dphidxsi[1] = 1.0;
  dphidxsi[2] = 0.0;
  dphideta[0] = -1.0;
  dphideta[1] = 0.0;
  dphideta[2] = 1.0;
}

void _e1c0_x(double *xsi) {
  xsi[0] = -1.0;
  xsi[1] = 1.0;
}

void _e1c0_phi(double xsi, double *phi) {
  phi[0] = (1 - xsi) / 2.0;
  phi[1] = (1 + xsi) / 2.0;
}

void _e1c0_dphidx(double xsi, double *dphidxsi) {
  dphidxsi[0] = -0.5;
  dphidxsi[1] = 0.5;
}

discrete *discreteCreate(int n, elementType type) {
  discrete *theSpace = malloc(sizeof(discrete));
  if (type == FEM_QUAD && n == 4) {
    theSpace->n = 4;
    theSpace->x2 = _q1c0_x;
    theSpace->phi2 = _q1c0_phi;
    theSpace->dphi2dx = _q1c0_dphidx;
  } else if (type == FEM_TRIANGLE && n == 3) {
    theSpace->n = 3;
    theSpace->x2 = _p1c0_x;
    theSpace->phi2 = _p1c0_phi;
    theSpace->dphi2dx = _p1c0_dphidx;
  } else if (type == FEM_EDGE && n == 2) {
    theSpace->n = 2;
    theSpace->x = _e1c0_x;
    theSpace->phi = _e1c0_phi;
    theSpace->dphidx = _e1c0_dphidx;
  } else
    Error("Cannot create such a discrete space !");
  return theSpace;
}

void discreteFree(discrete *theSpace) { free(theSpace); }

void discreteXsi(discrete *mySpace, double *xsi) { mySpace->x(xsi); }

void discretePhi(discrete *mySpace, double xsi, double *phi) { mySpace->phi(xsi, phi); }

void discreteXsi2(discrete *mySpace, double *xsi, double *eta) { mySpace->x2(xsi, eta); }

void discretePhi2(discrete *mySpace, double xsi, double eta, double *phi) { mySpace->phi2(xsi, eta, phi); }

void discreteDphi2(discrete *mySpace, double xsi, double eta, double *dphidxsi, double *dphideta) { mySpace->dphi2dx(xsi, eta, dphidxsi, dphideta); }

void discretePrint(discrete *mySpace) {
  int i, j;
  int n = mySpace->n;
  double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];

  discreteXsi2(mySpace, xsi, eta);
  for (i = 0; i < n; i++) {

    discretePhi2(mySpace, xsi[i], eta[i], phi);
    discreteDphi2(mySpace, xsi[i], eta[i], dphidxsi, dphideta);

    for (j = 0; j < n; j++) {
      printf("(xsi=%+.1f,eta=%+.1f) : ", xsi[i], eta[i]);
      printf(" phi(%d)=%+.1f", j, phi[j]);
      printf("   dphidxsi(%d)=%+.1f", j, dphidxsi[j]);
      printf("   dphideta(%d)=%+.1f \n", j, dphideta[j]);
    }
    printf(" \n");
  }
}
