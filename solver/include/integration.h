#pragma once

/**
 * @brief Structure containing the integration rule for a given element type (e.g. triangles or quads)
 * @param n Number of integration points
 * @param xsi Array of size n containing the xsi coordinates of the integration points
 * @param eta Array of size n containing the eta coordinates of the integration points
 * @param weight Array of size n containing the weights of the integration points
 */
typedef struct {
  int n;
  const double *xsi;
  const double *eta;
  const double *weight;
} integration;
