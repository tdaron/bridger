
#include <problem.h>
#include <solver.h>
#include <stdio.h>

void ss_init() { printf("[*] Solver Loaded\n"); }

double *compute_solution(char *filename, problem **prob, geo **geom, int nTankEdges, int* TankEdges, double tankWeight) {
  geo *theGeometry = geoMeshRead(filename);
  if (geom != NULL)
    *geom = theGeometry;
  // geoMeshPrint(theGeometry);

  double E = 211.e9;
  double nu = 0.3;
  double rho = 7.85e3;
  double g = 9.81;

  // Flag that determines if the matrix is to be banded or full
  int makeBanded = 1;

  problem *theProblem =
      elasticityCreate(theGeometry, E, nu, rho, g, PLANAR_STRAIN, makeBanded);
  if (prob != NULL)
    *prob = theProblem;
  // elasticityAddBoundaryCondition(theProblem,"Symmetry",DIRICHLET_X,0.0);
  // elasticityAddBoundaryCondition(theProblem,"Bottom",DIRICHLET_Y,0.0);
  // elasticityAddBoundaryCondition(theProblem,"Top",NEUMANN_Y,-1e4);

  // elasticityAddBoundaryCondition(theProblem,"PillarBottom_0",DIRICHLET_Y,0.0);
  // elasticityAddBoundaryCondition(theProblem,"PillarBottom_1",DIRICHLET_Y,0.0);
  // elasticityAddBoundaryCondition(theProblem,"PillarBottom_2",DIRICHLET_Y,0.0);
  // elasticityAddBoundaryCondition(theProblem,"PillarBottom_3",DIRICHLET_Y,0.0);

  // elasticityAddBoundaryCondition(theProblem,"PillarBottom_0",DIRICHLET_X,0.0);
  // elasticityAddBoundaryCondition(theProblem,"PillarBottom_1",DIRICHLET_X,0.0);
  // elasticityAddBoundaryCondition(theProblem,"PillarBottom_2",DIRICHLET_X,0.0);
  // elasticityAddBoundaryCondition(theProblem,"PillarBottom_3",DIRICHLET_X,0.0);

  // elasticityAddBoundaryCondition(theProblem,"LeftCorner",DIRICHLET_X,0.0);
  // elasticityAddBoundaryCondition(theProblem,"RightCorner",DIRICHLET_X,0.0);

  // elasticityAddBoundaryCondition(theProblem,"LeftCorner",DIRICHLET_Y,0.0);
  // elasticityAddBoundaryCondition(theProblem,"RightCorner",DIRICHLET_Y,0.0);

  elasticityAddBoundaryCondition(theProblem, "Pillar1", DIRICHLET_Y, 0.0);
  elasticityAddBoundaryCondition(theProblem, "Pillar2", DIRICHLET_Y, 0.0);
  elasticityAddBoundaryCondition(theProblem, "Pillar3", DIRICHLET_Y, 0.0);
  elasticityAddBoundaryCondition(theProblem, "Pillar4", DIRICHLET_Y, 0.0);

  elasticityAddBoundaryCondition(theProblem, "Pillar1", DIRICHLET_X, 0.0);
  elasticityAddBoundaryCondition(theProblem, "Pillar2", DIRICHLET_X, 0.0);
  elasticityAddBoundaryCondition(theProblem, "Pillar3", DIRICHLET_X, 0.0);
  elasticityAddBoundaryCondition(theProblem, "Pillar4", DIRICHLET_X, 0.0);

  elasticityAddBoundaryCondition(theProblem, "Extremity0", DIRICHLET_X, 0.0);
  elasticityAddBoundaryCondition(theProblem, "Extremity1", DIRICHLET_X, 0.0);

  elasticityAddBoundaryCondition(theProblem, "Extremity1", DIRICHLET_Y, 0.0);
  elasticityAddBoundaryCondition(theProblem, "Extremity0", DIRICHLET_Y, 0.0);

  
  elasticityAssembleNeumannExplicit(theProblem, TankEdges, nTankEdges, tankWeight);


  double *theSoluce = elasticitySolve(theProblem, makeBanded);

  // Reorder Solution
  double *solution_reorder =
      malloc(2 * theProblem->geometry->theNodes->nNodes * sizeof(double));
  for (int i = 0; i < 2 * theProblem->geometry->theNodes->nNodes; i++) {
    solution_reorder[2 * theProblem->renumNew2Old[i / 2] + i % 2] =
        theSoluce[i];
  }
  theSoluce = solution_reorder;
  nodes *theNodes = theGeometry->theNodes;
  double deformationFactor = 1e5;
  double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
  double *forcesX = malloc(theNodes->nNodes * sizeof(double));
  double *forcesY = malloc(theNodes->nNodes * sizeof(double));

  for (int i = 0; i < theNodes->nNodes; i++) {
    theNodes->X[i] += theSoluce[2 * i + 0] * deformationFactor;
    theNodes->Y[i] += theSoluce[2 * i + 1] * deformationFactor;
    normDisplacement[i] = sqrt(theSoluce[2 * i + 0] * theSoluce[2 * i + 0] +
                               theSoluce[2 * i + 1] * theSoluce[2 * i + 1]);
    // forcesX[i] = theForces[2*i+0];
    // forcesY[i] = theForces[2*i+1];
  }
  double hMin = vecMin(normDisplacement, theNodes->nNodes);
  double hMax = vecMax(normDisplacement, theNodes->nNodes);
  printf(" ==== Minimum displacement          : %14.7e [m] \n", hMin);
  printf(" ==== Maximum displacement          : %14.7e [m] \n", hMax);

  return theSoluce;
}
