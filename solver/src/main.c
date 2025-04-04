#include <problem.h>
#include <solver.h>
#include <stdio.h>

int main() {

  problem *theProblem;
  geo *theGeometry;
  double *theSoluce =
      compute_solution("../data/mesh.txt", &theProblem, &theGeometry, 0, NULL, 0);

  nodes *theNodes = theGeometry->theNodes;
  double deformationFactor = 1e5;
  double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));

  for (int i = 0; i < theNodes->nNodes; i++) {
    theNodes->X[i] += theSoluce[2 * i + 0] * deformationFactor;
    theNodes->Y[i] += theSoluce[2 * i + 1] * deformationFactor;
    normDisplacement[i] = sqrt(theSoluce[2 * i + 0] * theSoluce[2 * i + 0] +
                               theSoluce[2 * i + 1] * theSoluce[2 * i + 1]);
  }
}
