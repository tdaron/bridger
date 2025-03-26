#include <stdio.h>
#include <solver.h>

#include <problem.h>

double fun(double x, double y)
{
    return 1;
}

int main() {
    ss_init();
    printf("Hello, World\n");

    // geo* theGeometry = geoMeshRead("../data/elasticity.txt");
    geo* theGeometry = geoMeshRead("../data/mesh.txt");
    // geoMeshPrint(theGeometry);


    double E   = 211.e9;
    double nu  = 0.3;
    double rho = 7.85e3;
    double g   = 9.81;

    problem* theProblem = elasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRAIN);

    // elasticityAddBoundaryCondition(theProblem,"Symmetry",DIRICHLET_X,0.0);
    // elasticityAddBoundaryCondition(theProblem,"Bottom",DIRICHLET_Y,0.0);
    // elasticityAddBoundaryCondition(theProblem,"Top",NEUMANN_Y,-1e4);

    elasticityAddBoundaryCondition(theProblem,"Pillar1",DIRICHLET_Y,0.0);
    elasticityAddBoundaryCondition(theProblem,"Pillar2",DIRICHLET_Y,0.0);

    fprintf(stdout, "Debug 1\n");

    // elasticityPrint(theProblem);

    fprintf(stdout, "Debug 2\n");

    double *theSoluce = elasticitySolve(theProblem);

    fprintf(stdout, "Debug 3\n");
    double *theForces = elasticityForces(theProblem);

    fprintf(stdout, "Debug 4\n");
    double area = elasticityIntegrate(theProblem, fun);

    fprintf(stdout, "Debug 5\n");

    nodes *theNodes = theGeometry->theNodes;
    double deformationFactor = 1e5;
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    double *forcesX = malloc(theNodes->nNodes * sizeof(double));
    double *forcesY = malloc(theNodes->nNodes * sizeof(double));

    for (int i=0; i<theNodes->nNodes; i++){
        theNodes->X[i] += theSoluce[2*i+0]*deformationFactor;
        theNodes->Y[i] += theSoluce[2*i+1]*deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] +
                                    theSoluce[2*i+1]*theSoluce[2*i+1]);
        forcesX[i] = theForces[2*i+0];
        forcesY[i] = theForces[2*i+1]; }

    double hMin = vecMin(normDisplacement,theNodes->nNodes);
    double hMax = vecMax(normDisplacement,theNodes->nNodes);
    printf(" ==== Minimum displacement          : %14.7e [m] \n",hMin);
    printf(" ==== Maximum displacement          : %14.7e [m] \n",hMax);

    // double theGlobalForce[2] = {0, 0};
    // for (int i=0; i<theProblem->geometry->theNodes->nNodes; i++) {
    //     theGlobalForce[0] += theForces[2*i+0];
    //     theGlobalForce[1] += theForces[2*i+1]; }
    // printf(" ==== Global horizontal force       : %14.7e [N] \n",theGlobalForce[0]);
    // printf(" ==== Global vertical force         : %14.7e [N] \n",theGlobalForce[1]);
    // printf(" ==== Weight                        : %14.7e [N] \n", area * rho * g);

}
