#include "../include/problem.h"


int xycomp(const void *a, const void *b)
{
    // a, b are addresses of elements in the array "x[]",
    // each element x[i] has type "double *".
    double *valA = *(double**)a;   // Unpack the pointer
    double *valB = *(double**)b;   // Unpack the pointer

    // Now valA[0] and valB[0] are the actual coordinate
    if (valA[0] < valB[0]) return -1;
    if (valA[0] > valB[0]) return 1;
    return 0;
}

void problemXRenumber(problem *theProblem) {
    geo *problemGeomtery = theProblem->geometry;
    int i;
    double **x = malloc(sizeof(double *) * problemGeomtery->theNodes->nNodes);
    for (i = 0; i < problemGeomtery->theNodes->nNodes; i++) {
        x[i] = malloc(sizeof(double) * 2);
        x[i][0] = problemGeomtery->theNodes->X[i];
        x[i][1] = i;
    }
    qsort(x, problemGeomtery->theNodes->nNodes, sizeof(double *), xycomp);
    for (i = 0; i < problemGeomtery->theNodes->nNodes; i++) {
        theProblem->renumOld2New[(int)x[i][1]] = i;
        theProblem->renumNew2Old[i] = (int)x[i][1];
        // theProblem->renumOld2New[i] = i;
        // theProblem->renumNew2Old[i] = i;
        free(x[i]);
    }
    // Print the renumbering
    // printf("Renumbering : \n");
    // for (i = 0; i < problemGeomtery->theNodes->nNodes; i++) {
    //     fprintf(stderr, "%d ", theProblem->renumOld2New[i]);
    // }
    free(x);
}

int csr_comp(const void *a, const void *b)
{
    // a, b are addresses of elements in the array "x[]",
    // each element x[i] has type "double *".
    int **valA =  *(int***)a;   // Unpack the pointer
    int **valB = *(int***)b;   // Unpack the pointer

    if (valA[0][0] < valB[0][0]) return -1;
    if (valA[0][0] > valB[0][0]) return 1;
    if (valA[0][1] < valB[0][1]) return -1;
    return 1;
}

int** problemCOO(problem *theProblem) {
    int nNodes = theProblem->geometry->theNodes->nNodes;
    int nLocal = theProblem->geometry->theElements->nLocalNode;
    int* coo_adj_elems = malloc(2 * nNodes * nLocal * nLocal * sizeof(int));
    int** coo_adj = malloc(nNodes * nLocal * nLocal * sizeof(int*));
    for (int i = 0; i < nNodes; i++) {
        for (int j = 0; j < nLocal; j++) {
            for (int k = 0; k < nLocal; k++) {
                *coo_adj = coo_adj_elems;
                *coo_adj_elems = theProblem->geometry->theElements->elem[i * nLocal + j];
                *(coo_adj_elems + 1) = theProblem->geometry->theElements->elem[i * nLocal + k];
                coo_adj_elems += 2;
            }
        }
    }

    // Sort the COO matrix
    qsort(coo_adj, nNodes * nLocal * nLocal, sizeof(int*), csr_comp);

    // Now we have a list of pointers to tuples such that the
    // first element of the tuple is the row index and
    // the second element is the column index and where the
    // list is sorted by row index and then by column index.
    // Note : the numbers in memory have not been changed, only
    // the pointers to pointers to them have been sorted in
    // coo_adj.

    return coo_adj;
}

int** csr_from_coo(int **coo_adj, int nNodes, int nLocal) {
    int* col_ind = malloc(nNodes * sizeof(int));
    int* row_ptr = malloc((nNodes + 1) * sizeof(int));

    int prev_row = coo_adj[0][0];
    int* prev = NULL;
    for (int i = 0; i < nNodes * nLocal * nLocal; i++) {
        // Check for duplicates
        if (prev != NULL && prev[0] == coo_adj[i][0] && prev[1] == coo_adj[i][1]) {
            continue;
        }
        col_ind[i] = coo_adj[i][1];
        // Set all absent rows that were skipped to the
        // current column index.
        for (int j = prev_row; j < coo_adj[i][0]; j++) {
            row_ptr[j + 1] = i;
        }
        prev_row = coo_adj[i][0];
        prev = coo_adj[i];
    }

    // Make sure to update the last rows that might be empty
    for (int j = prev_row + 1 ; j < nNodes + 1; j++) {
        row_ptr[j] = nNodes + 1;
    }

    int** csr = malloc(2 * sizeof(int*));
    *csr = row_ptr;
    *(csr + 1) = col_ind;

    return csr;
}
