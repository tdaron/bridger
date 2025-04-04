#include "../include/problem.h"
#include <stdlib.h>
#include <stdbool.h>



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

int csr_comp(const void *a, const void *b) {
    // a and b are pointers to elements of coo_adj.
    // Each element of coo_adj is an int*.
    // So, we cast 'a' and 'b' to int** and dereference once to get the int* pointers.
    int *ptrA = *(int**)a;
    int *ptrB = *(int**)b;

    // Now ptrA and ptrB point to the pairs of integers (row, col) in coo_adj_elems
    int rowA = ptrA[0];
    int colA = ptrA[1];
    int rowB = ptrB[0];
    int colB = ptrB[1];

    // Compare rows first
    if (rowA < rowB) return -1;
    if (rowA > rowB) return 1;

    // Rows are equal, compare columns
    if (colA < colB) return -1;
    if (colA > colB) return 1; // Return 1 only if colA > colB

    // Rows and columns are equal
    return 0; // Important: Return 0 for equal elements
}

int* to_free_adj_elems = NULL;
int** to_free_adj = NULL;

int** problemCOO(problem *theProblem) {
    int nElems = theProblem->geometry->theElements->nElem;
    int nLocal = theProblem->geometry->theElements->nLocalNode;

    int num_pairs = nElems * nLocal * nLocal;

    int* coo_adj_elems = calloc(2 * num_pairs, sizeof(int));

    int** coo_adj = malloc(num_pairs * sizeof(int*));

    to_free_adj = coo_adj;
    to_free_adj_elems = coo_adj_elems;

    int* current_elem_ptr = coo_adj_elems;
    int current_adj_index = 0;

    for (int i = 0; i < nElems; i++) {
        for (int j = 0; j < nLocal; j++) {
            for (int k = 0; k < nLocal; k++) {
                // Assign the pointer in coo_adj[current_adj_index] to point
                // to the current location in coo_adj_elems
                coo_adj[current_adj_index] = current_elem_ptr;

                // Fill the data (row, col) in coo_adj_elems
                *current_elem_ptr = theProblem->geometry->theElements->elem[i * nLocal + j];       // Put row here
                *(current_elem_ptr + 1) = theProblem->geometry->theElements->elem[i * nLocal + k]; // Put col here

                current_elem_ptr += 2;
                current_adj_index++;
            }
        }
    }

    // Sort the array of pointers `coo_adj` based on the values of
    // the pairs in `coo_adj_elems` in increasing order.
    qsort(coo_adj, num_pairs, sizeof(int*), csr_comp);

    return coo_adj;
}



int** csr_from_coo(int **coo_adj, int nNodes, int nElems, int nLocal) {
    // int* col_ind = malloc(nNodes * sizeof(int));
    // int* row_ptr = malloc((nNodes + 1) * sizeof(int));
    int* col_ind = calloc(nElems * nLocal * nLocal, sizeof(int));
    int* row_ptr = calloc(nNodes + 1, sizeof(int));

    int prev_row = coo_adj[0][0];
    int* prev = NULL;
    int skipped = 0;
    for (int i = 0; i < nElems * nLocal * nLocal; i++) {
        // Check for duplicates
        if (prev != NULL && prev[0] == coo_adj[i][0] && prev[1] == coo_adj[i][1]) {
            skipped++;
            continue;
        }
        col_ind[i - skipped] = coo_adj[i][1];
        // Set all absent rows that were skipped to the
        // current column index.
        for (int j = prev_row; j < coo_adj[i][0]; j++) {
            row_ptr[j + 1] = i - skipped;}
        prev_row = coo_adj[i][0];
        prev = coo_adj[i];
    }
    // Make sure to update the last rows that might be empty
    for (int j = prev_row + 1 ; j < nNodes + 1; j++) {
        row_ptr[j] = nElems * nLocal * nLocal - skipped;
    }
    // Print size of col_ind and row_ptr and nNodes
    // fprintf(stderr, "Size of col_ind : %d\n", nElems * nLocal * nLocal - skipped);
    // fprintf(stderr, "Size of row_ptr : %d\n", nNodes + 1);
    // fprintf(stderr, "nNodes : %d\n", nNodes);

    int** csr = malloc(2 * sizeof(int*));
    *csr = row_ptr;
    *(csr + 1) = col_ind;

    return csr;
}

int* csr_row_glob = NULL;

int deg_sort(const void *a, const void *b)
{
    if (csr_row_glob == NULL) {
        fprintf(stderr, "Error : csr_row_glob is NULL in degree comparison function 🤡\n");
        return 0;
    }
    int idx1 = *(int*)a;
    int idx2 = *(int*)b;

    int deg1 = csr_row_glob[idx1 + 1] - csr_row_glob[idx1];
    int deg2 = csr_row_glob[idx2 + 1] - csr_row_glob[idx2];

    return deg1 - deg2;
}

// Implementation of the Reverse Cuthill-McKee algorithm for renumbering
// the nodes of a mesh to reduce the bandwidth of the matrix. Directly
// adapted from the original paper found at :
// https://dl.acm.org/doi/pdf/10.1145/800195.805928
// except for the order reversal at the end,
void rcm(problem* theProblem) {
    int nNodes = theProblem->geometry->theNodes->nNodes;
    int nElem = theProblem->geometry->theElements->nElem;
    int nLocal = theProblem->geometry->theElements->nLocalNode;

    int** coo_adj = problemCOO(theProblem);
    int* coo_adj_elems = coo_adj[0];

    int** csr = csr_from_coo(coo_adj, nNodes, nElem, nLocal);
    int* row_ptr = csr[0];
    int* col_ptr = csr[1];

    // Print the CSR col_ptr
    // for (int i = 0; i < nNodes; i++) {
    //     for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
    //         fprintf(stderr, "%d ", col_ptr[j]);
    //     }
    // }

    // Print the CSR row_ptr
    // for (int i = 0; i < nNodes + 1; i++) {
    //     fprintf(stderr, "%d ", row_ptr[i]);
    // }

    csr_row_glob = row_ptr;

    // Determine the node with the smallest degree
    int min_deg_idx = -1;
    int min_deg = -1;
    int deg;
    for (int i = 0; i < nNodes - 1; i++) {
        deg = row_ptr[i + 1] - row_ptr[i];
        if (min_deg == -1) {
            min_deg = deg;
            min_deg_idx = i;
        } else if (deg < min_deg) {
            min_deg = deg;
            min_deg_idx = i;
        }
    }

    // fprintf(stderr, "Min degree index : %d\n", min_deg_idx);
    // fprintf(stderr, "Min degree : %d\n", min_deg);

    int* queue = malloc(nNodes * sizeof(int));
    int* order = queue;
    // int* order = malloc(nNodes * sizeof(int));
    bool* visited = malloc(nNodes * sizeof(bool));
    for (int i = 0; i < nNodes; i++) {
        visited[i] = false;
    }

    queue[0] = min_deg_idx;
    // fprintf(stderr, "Queue : %d\n", queue[0]);
    // queue++;

    // Number of nodes that have been ordered,
    // used to known when to stop the BFS.
    int nVisited = 0;

    // Number of nodes in the queue
    int nInQueue = 0;

    // Number of nodes in the queue in the previous iteration,
    // used to know which part of the queue is new and so which
    // part of the queue to sort based on degree.
    int old_nInQueue = 0;

    int skipped = 0;
    visited[min_deg_idx] = true;
    // fprintf(stderr, "Visited : %d\n", min_deg_idx);
    for (int i = row_ptr[min_deg_idx]; i < row_ptr[min_deg_idx + 1]; i++) {
        // fprintf(stderr, "Col : %d\n", col_ptr[i]);
        // if (!visited[col_ptr[i]]) {
            // fprintf(stderr, "Col2 : %d\n", col_ptr[i]);
            // fprintf(stderr, "nInQueue : %d\n", nInQueue);
        queue[nInQueue] = col_ptr[i];
        visited[col_ptr[i]] = true;
        nInQueue++;
        // if (col_ptr[i] == 1726) {
            // fprintf(stderr, "Col : %d\n", col_ptr[i]);
        // }
        // }
    }
    // fprintf(stderr, "Queue : %d\n", queue[0]);
    qsort(queue, nInQueue, sizeof(int), deg_sort);
    // Print the queue
    // for (int i = 0; i < nInQueue; i++) {
    //     fprintf(stderr, "%d ", queue[i]);
    // }
    while (nVisited < nNodes) {
        int node = queue[0];
        queue++;
        nInQueue--;
        old_nInQueue = nInQueue;
        // fprintf(stderr, "%d ", node);
        // if (node == 0) {
        //     exit(-1);
        // }

        // if (visited[node]) {
        //     continue;
        // }
        nVisited++;
        // Save the number of nodes that were in the queue,
        // in the previous iteration.
        // Iterate over the neighbors of the current node
        for (int i = row_ptr[node]; i < row_ptr[node + 1]; i++) {
            if (!visited[col_ptr[i]]) {
                // old_nInQueue is used here to determine the offset
                // from which we can start adding new nodes to the queue
                // without overwriting the old ones.
                // queue[old_nInQueue + i - row_ptr[node]] = col_ptr[i];
                queue[nInQueue] = col_ptr[i];
                // if (col_ptr[i] == 1726) {
                //     fprintf(stderr, "Col : %d\n", col_ptr[i]);
                // }
                visited[col_ptr[i]] = true;
                nInQueue++;
            }
        }
        // Sort the new nodes in the queue based on degree
        qsort(queue + old_nInQueue, nInQueue - old_nInQueue, sizeof(int), deg_sort);
    }

    // fprintf(stderr, "\n");
    // fprintf(stderr, "\n");
    // fprintf(stderr, "\n");
    // fprintf(stderr, "\n");
    for (int i = 0; i < nNodes; i++) {
        // Standard Cuthill-McKee
        // theProblem->renumOld2New[order[i]] = i;
        // theProblem->renumNew2Old[i] = order[i];

        // Reverse Cuthill-McKee
        theProblem->renumOld2New[order[nNodes - 1 - i]] = i;
        // fprintf(stderr, "%i ", order[nNodes - 1 - i]);
        theProblem->renumNew2Old[i] = order[nNodes - 1 - i];
    }

    // fprintf(stderr, "\n");

    // Print the renumbering
    // fprintf(stderr, "Renumbering : \n");
    // for (int i = 0; i < nNodes; i++) {
    //     fprintf(stderr, "%d ", theProblem->renumOld2New[i]);
    // }


    free(order);
    free(visited);
    free(row_ptr);
    free(col_ptr);
    free(csr);
    free(to_free_adj_elems);
    free(to_free_adj);
}
