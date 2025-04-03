#include <problem.h>
#include <solver.h>
#include <stdio.h>

int main() {

  problem *theProblem;
  geo *theGeometry;
  double *theSoluce =
      compute_solution("../data/mesh.txt", &theProblem, &theGeometry);

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

  FILE *file = fopen("solution.txt", "w");

  for (int i = 0; i < theProblem->geometry->theNodes->nNodes; i++) {

    // 10 digit precsion

    fprintf(file, "%14.7e %14.7e %d\n", theSoluce[2 * i + 0],
            theSoluce[2 * i + 1], i);
  }

  fclose(file);
  // double theGlobalForce[2] = {0, 0};
  // for (int i=0; i<theProblem->geometry->theNodes->nNodes; i++) {
  //     theGlobalForce[0] += theForces[2*i+0];
  //     theGlobalForce[1] += theForces[2*i+1]; }
  // printf(" ==== Global horizontal force       : %14.7e [N]
  // \n",theGlobalForce[0]); printf(" ==== Global vertical force         :
  // %14.7e [N] \n",theGlobalForce[1]); printf(" ==== Weight : %14.7e [N] \n",
  // area * rho * g);
}
