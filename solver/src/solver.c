#include "../include/problem.h"

double *elasticitySolve(problem *theProblem, int makeBanded) {

    fullSystemInit(theProblem->system);

    if (makeBanded) {
        elasticityAssembleElementsBand(theProblem);
    } else {
        elasticityAssembleElements(theProblem);
    }
    elasticityAssembleNeumann(theProblem);

    int *theConstrainedNodes = theProblem->constrainedNodes;
    for (int i = 0; i < theProblem->system->size; i++) {
      if (theConstrainedNodes[i] != -1) {
        double value = theProblem->conditions[theConstrainedNodes[i]]->value;
        if (makeBanded) {
          bandSystemConstrain(theProblem->system, i, value);
        } else {
          fullSystemConstrain(theProblem->system, i, value);
        }
      }
    }

    if (makeBanded) {
        bandSystemEliminate(theProblem->system);
    } else {
        fullSystemEliminate(theProblem->system);
    }

    for (int i = 0; i < theProblem->system->size; i++) {
        theProblem->soluce[i] = theProblem->system->B[i];
    }

    return theProblem->soluce;
}

double *fullSystemEliminate(linearSystem *mySystem) {
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

double *bandSystemEliminate(linearSystem *myBand) {
  double **A, *B, factor;
  int i, j, k, jend, size, band;
  A = myBand->A;
  B = myBand->B;
  size = myBand->size;
  band = myBand->band;

  // A completer :-)

  for (k = 0; k < size; k++) {
    if (fabs(A[k][k]) <= 1e-8) {
      printf("Pivot index %d  ", k);
      printf("Pivot value %e  ", A[k][k]);
      Error("Cannot eliminate with such a pivot");
    }
    jend = fmin(size, k + band);
    for (i = k + 1; i < jend; i++) {
      // Here we have swapped the usual A[i][k] for A[k][i]
      // because gaussian elimination preserves symmetry
      factor = A[k][i] / A[k][k];
      for (j = i; j < jend; j++) {
        A[i][j] -= A[k][j] * factor;
      }
      B[i] -= B[k] * factor;
    }
  }
  for (i = size - 1; i >= 0; i--) {
    factor = 0;
    jend = fmin(size, i + band);
    for (j = i + 1; j < jend; j++) {
      factor += A[i][j] * B[j];
    }
    B[i] = (B[i] - factor) / A[i][i];
  }

  return (myBand->B);
}
