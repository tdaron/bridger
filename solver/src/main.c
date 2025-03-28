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

    // Flag that determines if the matrix is to be banded or full
    int makeBanded = 0;

    problem* theProblem = elasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRAIN, makeBanded);

    // elasticityAddBoundaryCondition(theProblem,"Symmetry",DIRICHLET_X,0.0);
    // elasticityAddBoundaryCondition(theProblem,"Bottom",DIRICHLET_Y,0.0);
    // elasticityAddBoundaryCondition(theProblem,"Top",NEUMANN_Y,-1e4);

    elasticityAddBoundaryCondition(theProblem,"Pillar1",DIRICHLET_Y,0.0);
    elasticityAddBoundaryCondition(theProblem,"Pillar2",DIRICHLET_Y,0.0);
    elasticityAddBoundaryCondition(theProblem,"Pillar3",DIRICHLET_Y,0.0);
    elasticityAddBoundaryCondition(theProblem,"Pillar4",DIRICHLET_Y,0.0);

    elasticityAddBoundaryCondition(theProblem,"Left",DIRICHLET_X,0.0);
    elasticityAddBoundaryCondition(theProblem,"Right",DIRICHLET_X,0.0);

    double *theSoluce = elasticitySolve(theProblem, makeBanded);

    // Reorder Solution
    double * solution_reorder = malloc(2*theProblem->geometry->theNodes->nNodes*sizeof(double));
    for (int i=0; i<2*theProblem->geometry->theNodes->nNodes; i++) {
        solution_reorder[2 * theProblem->renumNew2Old[i / 2] + i % 2] = theSoluce[i];
    }
    theSoluce = solution_reorder;

    // Wrtie the solution to a file
    // FILE *file = fopen("solution.txt", "w");
    // for (int i=0; i<theProblem->geometry->theNodes->nNodes; i++) {
    //     // 10 digit precsion
    //     fprintf(file, "%14.7e %14.7e\n", theSoluce[2*i+0], theSoluce[2*i+1]);
    // }
    // fclose(file);

    // double *theForces = elasticityForces(theProblem);
    // double area = elasticityIntegrate(theProblem, fun);

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
        // forcesX[i] = theForces[2*i+0];
        // forcesY[i] = theForces[2*i+1];
    }

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
