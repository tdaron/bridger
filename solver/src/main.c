#include <problem.h>
#include <solver.h>
#include <stdio.h>

int main() {

  fprintf(stderr, "Hello World\n");
  problem *theProblem;
  geo *theGeometry;
  double *theSoluce =
      compute_solution("./data/mesh.txt", &theProblem, &theGeometry, 0, NULL, 0);

  fprintf(stderr, " ==== Problem solved \n");
  nodes *theNodes = theGeometry->theNodes;
  double deformationFactor = 1e5;
  int nNodes = theNodes->nNodes;
  double *normDisplacement = malloc(nNodes * sizeof(double));

  for (int i = 0; i < theNodes->nNodes; i++) {
    theNodes->X[i] += theSoluce[2 * i + 0] * deformationFactor;
    theNodes->Y[i] += theSoluce[2 * i + 1] * deformationFactor;
    normDisplacement[i] = sqrt(theSoluce[2 * i + 0] * theSoluce[2 * i + 0] +
                               theSoluce[2 * i + 1] * theSoluce[2 * i + 1]);
  }

  double hMin = vecMin(normDisplacement,nNodes);
  double hMax = vecMax(normDisplacement,nNodes);
  printf(" ==== Minimum displacement          : %14.7e [m] \n",hMin);
  printf(" ==== Maximum displacement          : %14.7e [m] \n",hMax);

  geoFree(theGeometry);
  // elasticityFree(theProblem);
  free(theGeometry);
  free(theSoluce);
}
