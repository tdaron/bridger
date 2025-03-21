#include "../include/problem.h"

double *elasticitySolve(problem *theProblem){

    //
    // A completer :-)
    //

    fullSystemInit(theProblem->system);

    elasticityAssembleElements(theProblem);
    // fprintf(stdout, "Assemble elements done\n");
    elasticityAssembleNeumann(theProblem);
    // fprintf(stdout, "Assemble Neumann done\n");


    int *theConstrainedNodes = theProblem->constrainedNodes;
    for (int i = 0; i < theProblem->system->size; i++) {
      if (theConstrainedNodes[i] != -1) {
        // fprintf(stdout, "Constraining node %d\n", i);
        double value = theProblem->conditions[theConstrainedNodes[i]]->value;
        fullSystemConstrain(theProblem->system, i, value);
      }
    }


    fullSystemEliminate(theProblem->system);

    for (int i = 0; i < theProblem->system->size; i++) {
        theProblem->soluce[i] = theProblem->system->B[i];
    }

    fprintf(stdout, "Done\n");

    //

     return theProblem->soluce;
}

double *fullSystemEliminate(fullSystem *mySystem) {
  double **A, *B, factor;
  int i, j, k, size;

  A = mySystem->A;
  B = mySystem->B;
  size = mySystem->size;

  /* Gauss elimination */

  for (k = 0; k < size; k++) {
    if (fabs(A[k][k]) <= 1e-16) {
      printf("Pivot index %d  ", k);
      printf("Pivot value %e  ", A[k][k]);
      Error("Cannot eliminate with such a pivot");
    }
    for (i = k + 1; i < size; i++) {
      factor = A[i][k] / A[k][k];
      for (j = k + 1; j < size; j++)
        A[i][j] = A[i][j] - A[k][j] * factor;
      B[i] = B[i] - B[k] * factor;
    }
  }

  /* Back-substitution */

  for (i = size - 1; i >= 0; i--) {
    factor = 0;
    for (j = i + 1; j < size; j++)
      factor += A[i][j] * B[j];
    B[i] = (B[i] - factor) / A[i][i];
  }

  return (mySystem->B);
}


double * elasticityForces(problem *theProblem)
{
    // // Initialize system first so we don't lose the matrix after assembly
    // fullSystemInit(theProblem->system);
    double *soluce = theProblem->soluce;

    // Copy solution
    double* savedSoluce = malloc(theProblem->system->size * sizeof(double));
    for (int i = 0; i < theProblem->system->size; i++) {
        savedSoluce[i] = soluce[i];
    }

    fullSystemInit(theProblem->system);

    // Reassemble the elements (matrix A and vector B will be populated here)
    elasticityAssembleElements(theProblem);
    // fprintf(stdout, "Assemble elements done\n");
    elasticityAssembleNeumann(theProblem);
    // fprintf(stdout, "Assemble Neumann done\n");

    // Now compute the residuals based on DOF indexing
    int size = theProblem->system->size;  // number_of_nodes
    for (int i = 0; i < size; i++) {
        // Make sure you reset the residual to 0 before summation
        theProblem->residuals[i]   = 0.0;

        for (int j = 0; j < size; j++) {
            theProblem->residuals[i]   += theProblem->system->A[i][j]
                                            * savedSoluce[j];
        }

        theProblem->residuals[i]   -= theProblem->system->B[i];
    }

    // Print A

    free(savedSoluce);


    return theProblem->residuals;
}
