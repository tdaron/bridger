#include "../include/problem.h"


void elasticityAssembleElements(problem *theProblem){
    linearSystem  *theSystem = theProblem->system;
    integration *theRule = theProblem->rule;
    discrete    *theSpace = theProblem->space;
    geo         *theGeometry = theProblem->geometry;
    nodes       *theNodes = theGeometry->theNodes;
    mesh        *theMesh = theGeometry->theElements;
    mesh        *theEdges = theGeometry->theEdges;
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;


    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];
            map[j] = theProblem->renumOld2New[map[j]];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
        }

        for (iInteg=0; iInteg < theRule->n; iInteg++) {
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            discretePhi2(theSpace,xsi,eta,phi);
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

            for (i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }
            for (i = 0; i < theSpace->n; i++) {
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] +
                                            dphidy[i] * c * dphidy[j]) * jac * weight;
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] +
                                            dphidy[i] * c * dphidx[j]) * jac * weight;
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] +
                                            dphidx[i] * c * dphidy[j]) * jac * weight;
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] +
                                            dphidx[i] * c * dphidx[j]) * jac * weight; }}
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; }}}
}

 // Calcule et assemble toutes les intégrales de ligne.
 // Ce sont les intégrales associées aux conditions de
 // Neumann de la formulation discrète.
void elasticityAssembleNeumann(problem *theProblem){
    linearSystem  *theSystem = theProblem->system;
    integration *theRule = theProblem->ruleEdge;
    discrete    *theSpace = theProblem->spaceEdge;
    geo         *theGeometry = theProblem->geometry;
    nodes       *theNodes = theGeometry->theNodes;
    mesh        *theEdges = theGeometry->theEdges;
    double x[2],y[2],phi[2], dphidx[2];
    int iBnd,iElem,iInteg,iEdge,i,j,d,map[2],mapU[2];
    int nLocal = 2;
    double *B  = theSystem->B;

    for(iBnd=0; iBnd < theProblem->nBoundaryConditions; iBnd++){
        boundaryCondition *theCondition = theProblem->conditions[iBnd];
        boundaryType type = theCondition->type;
        double value = theCondition->value;

        if (type == DIRICHLET_X || type == DIRICHLET_Y) {
            continue;
        }

        int shift = (type == NEUMANN_X) ? 0 : 1;

        // Iterate over the number of edges in the domain
        for (iEdge = 0; iEdge < theCondition->domain->nElem; iEdge++) {
            // Get the actual element index
            iElem = theCondition->domain->elem[iEdge];

            // For each edge, we only have the shape functions
            // associated with the two nodes of the edge. that
            // are non-zero.
            for (j = 0; j < nLocal; j++) {
                map[j] = theEdges->elem[iElem * nLocal + j];
                x[j] = theNodes->X[map[j]];
                y[j] = theNodes->Y[map[j]];
                map[j] = theProblem->renumOld2New[map[j]];
            }

            // Apply the integration rule to the edge.
            for (iInteg = 0; iInteg < theRule->n; iInteg++) {
                double xsi = theRule->xsi[iInteg];
                double weight = theRule->weight[iInteg];

                theSpace->phi(xsi, phi);
                theSpace->dphidx(xsi, dphidx);

                double dxdxsi = 0;
                for (i = 0; i < theSpace->n; i++) {
                    dxdxsi += x[i] * dphidx[i];
                }
                double jac = fabs(dxdxsi);
                for (i = 0; i < theSpace->n; i++) {
                    B[2 * map[i] + shift] += phi[i] * value * weight * jac;
                }
            }
        }
    }
}


void elasticityAssembleElementsBand(problem *theProblem){
    linearSystem  *theSystem = theProblem->system;
    integration *theRule = theProblem->rule;
    discrete    *theSpace = theProblem->space;
    geo         *theGeometry = theProblem->geometry;
    nodes       *theNodes = theGeometry->theNodes;
    mesh        *theMesh = theGeometry->theElements;
    mesh        *theEdges = theGeometry->theEdges;
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;


    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];
            map[j] = theProblem->renumOld2New[map[j]];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
        }

        for (iInteg=0; iInteg < theRule->n; iInteg++) {
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            discretePhi2(theSpace,xsi,eta,phi);
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

            for (i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }
            for (i = 0; i < theSpace->n; i++) {
                for(j = 0; j < theSpace->n; j++) {
                    if (mapX[j] <= mapX[i]) {
                        A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] +
                                                dphidy[i] * c * dphidy[j]) * jac * weight;
                    }
                    if (mapY[j] <= mapX[i]) {
                        A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] +
                                                dphidy[i] * c * dphidx[j]) * jac * weight;
                    }
                    if (mapX[j] <= mapY[i]) {
                        A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] +
                                                dphidx[i] * c * dphidy[j]) * jac * weight;
                    }
                    if (mapY[j] <= mapY[i]) {
                        A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] +
                                                dphidx[i] * c * dphidx[j]) * jac * weight;
                    }
                }
            }
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; }}}
}
