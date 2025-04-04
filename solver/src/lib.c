
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

  // elasticityAddBoundaryCondition(theProblem, "Tank", DIRICHLET_Y, 0.0);
  // elasticityAddBoundaryCondition(theProblem, "Tank", DIRICHLET_Y, 0.0);


  // Print max value of B vector
  // fprintf(stdout, "Max B value: %e\n", vecMax(theProblem->system->B, theProblem->system->size));


  // elasticityAssembleNeumannExplicit(theProblem, TankEdges, nTankEdges, tankWeight);


  double *theSoluce = elasticitySolve(theProblem, makeBanded, TankEdges, nTankEdges, tankWeight);

  // Reorder Solution
  double *solution_reorder =
      malloc(2 * theProblem->geometry->theNodes->nNodes * sizeof(double));
  for (int i = 0; i < 2 * theProblem->geometry->theNodes->nNodes; i++) {
    solution_reorder[2 * theProblem->renumNew2Old[i / 2] + i % 2] =
        theSoluce[i];
  }
  theSoluce = solution_reorder;



  return theSoluce;
}
