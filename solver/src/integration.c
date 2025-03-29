#include "../include/problem.h"

double elasticityIntegrate(problem *theProblem, double (*f)(double x, double y)){
    integration *theRule = theProblem->rule;
    geo         *theGeometry = theProblem->geometry;
    nodes       *theNodes = theGeometry->theNodes;
    mesh        *theMesh = theGeometry->theElements;
    discrete    *theSpace = theProblem->space;

    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,i,map[4];
    int nLocal = theMesh->nLocalNode;
    double value = 0.0;
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (i=0; i < nLocal; i++) {
            map[i]  = theMesh->elem[iElem*nLocal+i];
            x[i]    = theNodes->X[map[i]];
            y[i]    = theNodes->Y[map[i]];}
        for (iInteg=0; iInteg < theRule->n; iInteg++) {
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            discretePhi2(theProblem->space,xsi,eta,phi);
            discreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0;
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {
                dxdxsi += x[i]*dphidxsi[i];
                dxdeta += x[i]*dphideta[i];
                dydxsi += y[i]*dphidxsi[i];
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theProblem->space->n; i++) {
                value += phi[i] * f(x[i],y[i]) * jac * weight; }}}
    return value;
}

static const double _gaussQuad4Xsi[4] = {-0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Eta[4] = {0.577350269189626, -0.577350269189626, -0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Weight[4] = {1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000};
static const double _gaussTri3Xsi[3] = {0.166666666666667, 0.666666666666667, 0.166666666666667};
static const double _gaussTri3Eta[3] = {0.166666666666667, 0.166666666666667, 0.666666666666667};
static const double _gaussTri3Weight[3] = {0.166666666666667, 0.166666666666667, 0.166666666666667};
static const double _gaussEdge2Xsi[2] = {0.577350269189626, -0.577350269189626};
static const double _gaussEdge2Weight[2] = {1.000000000000000, 1.000000000000000};

integration *integrationCreate(int n, elementType type) {
  integration *theRule = malloc(sizeof(integration));
  if (type == FEM_QUAD && n == 4) {
    theRule->n = 4;
    theRule->xsi = _gaussQuad4Xsi;
    theRule->eta = _gaussQuad4Eta;
    theRule->weight = _gaussQuad4Weight;
  } else if (type == FEM_TRIANGLE && n == 3) {
    theRule->n = 3;
    theRule->xsi = _gaussTri3Xsi;
    theRule->eta = _gaussTri3Eta;
    theRule->weight = _gaussTri3Weight;
  } else if (type == FEM_EDGE && n == 2) {
    theRule->n = 2;
    theRule->xsi = _gaussEdge2Xsi;
    theRule->eta = NULL;
    theRule->weight = _gaussEdge2Weight;
  } else
    Error("Cannot create such an integration rule !");
  return theRule;
}

void integrationFree(integration *theRule) { free(theRule); }
