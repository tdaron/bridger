#pragma once
#include <problem.h>
void ss_init();
double *compute_solution(char *filename, problem **prob, geo **geom, int nTankEdges, int* TankEdges, double tankWeight);
double *compute_solution2(char *filename, problem **prob, geo **geom, int nTankEdges, int* TankEdges, double tankWeight);
